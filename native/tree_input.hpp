#pragma once

#include "alisim_bridge.h"
#include "tree/iqtree.h"
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace pf {
struct TreeInput {
    uint32_t tips = 0;
    bool rooted = false;
    std::vector<std::vector<std::pair<uint32_t, double>>> neighbors;
};
inline TreeInput& tree_input() { static TreeInput input; return input; }

// The Rust wrapper serializes calls. Only plain C values cross the ABI; native
// nodes are created through the actual tree's virtual factory and owned by it.
struct TreeSession {
    TreeSession(const pf_tree_edge* edges, size_t count, uint32_t tips) {
        tree_input() = TreeInput{};
        if (!edges) {
            if (count || tips) throw std::invalid_argument("Missing tree edges");
            return;
        }
        if (tips < 3 || tips > 10000 || (count != 2 * size_t(tips) - 3 && count != 2 * size_t(tips) - 2))
            throw std::invalid_argument("Invalid binary tree dimensions");
        TreeInput input;
        input.tips = tips;
        input.rooted = count == 2 * size_t(tips) - 2;
        input.neighbors.resize(count + 1);
        for (size_t i = 0; i < count; ++i) {
            const auto& edge = edges[i];
            if (edge.a >= input.neighbors.size() || edge.b >= input.neighbors.size() ||
                edge.a == edge.b || !std::isfinite(edge.length) || edge.length < 0)
                throw std::invalid_argument("Invalid tree edge");
            input.neighbors[edge.a].push_back({edge.b, edge.length});
            input.neighbors[edge.b].push_back({edge.a, edge.length});
        }
        for (size_t i = 0; i < input.neighbors.size(); ++i)
            if (input.neighbors[i].size() != (i < tips ? 1U : (input.rooted && i == tips ? 2U : 3U)))
                throw std::invalid_argument("Invalid tree node degree");
        // Validate connectivity before any recursive IQ-TREE operation.
        std::vector<bool> seen(input.neighbors.size(), false);
        std::vector<uint32_t> stack{tips};
        size_t visited = 0;
        while (!stack.empty()) {
            const auto node = stack.back();
            stack.pop_back();
            if (seen[node]) continue;
            seen[node] = true;
            ++visited;
            for (const auto& neighbor : input.neighbors[node])
                if (!seen[neighbor.first]) stack.push_back(neighbor.first);
        }
        if (visited != input.neighbors.size())
            throw std::invalid_argument("Tree is disconnected or contains a cycle");
        tree_input() = std::move(input);
    }
    ~TreeSession() { tree_input() = TreeInput{}; }
    TreeSession(const TreeSession&) = delete;
    TreeSession& operator=(const TreeSession&) = delete;
};

inline void read_tree(IQTree& tree, const char* filename, bool& is_rooted) {
    const auto& input = tree_input();
    if (!input.tips) {
        tree.readTree(filename, is_rooted);
        return;
    }
    if (tree.root || tree.isSuperTree())
        throw std::invalid_argument("Direct tree input requires a fresh, unpartitioned tree");

    std::vector<std::unique_ptr<Node>> nodes;
    nodes.reserve(input.neighbors.size());
    for (size_t i = 0; i < input.neighbors.size(); ++i) {
        const auto name = i < input.tips ? "T" + std::to_string(i + 1) : "";
        nodes.emplace_back(tree.newNode(-1, name.c_str()));
    }

    // Match Tree::newick's traversal and MTree::parseFile's postorder edge
    // attachment, leaf IDs and root selection. Neighbor order affects RNG use.
    struct Visit { uint32_t node; uint32_t parent; double length; bool exit; };
    const uint32_t no_parent = UINT32_MAX;
    std::vector<Visit> stack{{input.tips, no_parent, 0, false}};
    Node* root = nullptr;
    int next_leaf = 0;
    while (!stack.empty()) {
        const auto visit = stack.back();
        stack.pop_back();
        auto* node = nodes[visit.node].get();
        if (visit.exit) {
            if (visit.parent != no_parent) {
                auto* parent = nodes[visit.parent].get();
                parent->addNeighbor(node, visit.length);
                node->addNeighbor(parent, visit.length);
            }
            continue;
        }
        if (visit.node < input.tips) {
            node->id = next_leaf++;
            if (!root) root = node;
        }
        stack.push_back({visit.node, visit.parent, visit.length, true});
        const auto& neighbors = input.neighbors[visit.node];
        for (auto it = neighbors.rbegin(); it != neighbors.rend(); ++it)
            if (it->first != visit.parent)
                stack.push_back({it->first, visit.node, it->second, false});
    }
    if (!input.rooted) for (const auto& neighbor : input.neighbors[input.tips])
        if (neighbor.first < input.tips) {
            root = nodes[neighbor.first].get();
            break;
        }

    std::unique_ptr<Node> root_leaf;
    if (input.rooted) {
        // Match MTree's explicit virtual root leaf for a bifurcating Newick root.
        root_leaf.reset(tree.newNode(input.tips, ROOT_NAME));
        auto* crown = nodes[input.tips].get();
        crown->addNeighbor(root_leaf.get(), 0.0);
        root_leaf->addNeighbor(crown, 0.0);
        root = root_leaf.get();
    }
    tree.root = root;
    tree.leafNum = input.tips + (input.rooted ? 1 : 0);
    tree.rooted = is_rooted = input.rooted;
    root_leaf.release();
    for (auto& node : nodes) node.release();
    tree.initializeTree();
    if (input.rooted) tree.computeBranchDirection();
}
} // namespace pf
