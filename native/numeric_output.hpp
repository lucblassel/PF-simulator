#pragma once

#include "main/alisim.h"
#include <memory>
#include <algorithm>

// Numeric transport used by the simulator; FASTA is produced only after acceptance.
namespace pf {
struct RawRow {
    std::string name;
    std::vector<uint8_t> states;
};
struct OriginalRow {
    std::string name;
    std::vector<short int> states;
};
struct RawState {
    bool active = false;
    std::vector<RawRow> rows;
    std::vector<OriginalRow> originals;
    std::vector<uint8_t> alphabet;
    std::vector<uint8_t> canonical_codes;
};
inline RawState& raw_state() { static RawState state; return state; }
inline bool raw_enabled() { return raw_state().active; }
struct RawSession {
    explicit RawSession(bool active) { raw_state() = RawState{}; raw_state().active = active; }
    ~RawSession() { raw_state() = RawState{}; }
};

inline void capture_raw(const std::string& name, const std::vector<short int>& states,
                        int length, const std::vector<std::string>& mapping) {
    auto& data = raw_state();
    if (data.alphabet.empty()) {
        for (const auto& symbol : mapping) {
            if (symbol.size() != 1) throw std::runtime_error("Numeric output supports DNA/protein only");
            data.alphabet.push_back(static_cast<uint8_t>(symbol[0]));
        }
    }
    if (data.canonical_codes.empty()) {
        if (data.alphabet.size() > 256) throw std::runtime_error("Numeric alphabet exceeds one byte");
        for (const auto symbol : data.alphabet)
            data.canonical_codes.push_back(static_cast<uint8_t>(std::find(data.alphabet.begin(), data.alphabet.end(), symbol) - data.alphabet.begin()));
    }
    RawRow row{name, {}};
    row.states.reserve(length);
    for (int i = 0; i < length; ++i) {
        const auto state = states.at(i);
        if (state < 0 || state > 255 || state >= mapping.size()) {
            throw std::runtime_error("Numeric output state does not fit its byte alphabet");
        }
        // Canonicalize aliases so byte equality matches exported character equality.
        row.states.push_back(data.canonical_codes[state]);
    }
    data.rows.push_back(std::move(row));
}

inline std::string sequence_name(const Node* node) {
    return node->name.empty() ? std::to_string(node->id) : node->name;
}

inline void capture_retained_rows(Node* node, Node* dad, bool skip_leaves,
                                  int length, const std::vector<std::string>& mapping) {
    const auto& params = Params::getInstance();
    if (!(node->isLeaf() && skip_leaves) &&
        ((node->isLeaf() && node->name != ROOT_NAME) ||
         (params.alisim_write_internal_sequences && params.alisim_insertion_ratio + params.alisim_deletion_ratio > 0)))
        capture_raw(sequence_name(node), node->sequence->sequence_chunks[0], length, mapping);
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        capture_retained_rows((*it)->node, node, skip_leaves, length, mapping);
    }
}

// Use AliSim's genome reconstruction and preserve the original traversal and
// insertion-history updates, but retain numeric rows instead of text states.
inline void capture_indel_rows(AliSimulator* simulator, int length,
                               const std::vector<std::string>& mapping) {
    std::unique_ptr<GenomeTree> genome;
    Insertion* previous = nullptr;
    const int step = simulator->params->rebuild_indel_history_param * simulator->tree->leafNum;
    int threshold = step;
    int line = 0;
    for (auto& original : raw_state().originals) {
        ++line;
        auto* node = simulator->map_seqname_node.at(original.name);
        const int size = static_cast<int>(original.states.size());
        if (!genome) {
            genome = std::make_unique<GenomeTree>();
            genome->buildGenomeTree(node->sequence->insertion_pos, size, true);
        } else if (node->sequence->insertion_pos->next) {
            if (line >= threshold) {
                for (Insertion* insertion = node->sequence->insertion_pos; insertion; insertion = insertion->next)
                    insertion->genome_nodes.clear();
                genome = std::make_unique<GenomeTree>();
                genome->buildGenomeTree(node->sequence->insertion_pos, size, true);
                threshold += step;
            } else {
                genome->updateGenomeTree(previous, node->sequence->insertion_pos);
            }
        } else {
            genome = std::make_unique<GenomeTree>(size);
        }
        previous = node->sequence->insertion_pos;
        node->sequence->insertion_pos = nullptr;
        const auto aligned = genome->exportNewGenome(original.states, length, simulator->tree->aln->STATE_UNKNOWN);
        capture_raw(original.name, aligned, length, mapping);
    }
}
}  // namespace pf
