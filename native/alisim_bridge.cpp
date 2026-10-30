#include "alisim_bridge.h"
#include "memory_io.hpp"
#include "tree_input.hpp"
#include "numeric_output.hpp"
#include "main/alisim.h"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace {

// Redirect native informational streams without touching stderr's file
// descriptor: Rust's progress bar continues drawing while AliSim is running.
// Errors on std::cerr remain visible, including upstream exit() error paths.
struct OutputRedirect {
    std::ofstream output;
    std::streambuf* previous_cout;
    std::streambuf* previous_clog;

    explicit OutputRedirect(const char* log_path)
        : output(log_path ? log_path : "/dev/null", std::ios::app) {
        if (!output) throw std::runtime_error("Cannot open AliSim log output");
        previous_cout = std::cout.rdbuf(output.rdbuf());
        previous_clog = std::clog.rdbuf(output.rdbuf());
    }

    void finish() {
        output.flush();
        if (!output) throw std::runtime_error("Cannot write AliSim log output");
    }

    ~OutputRedirect() {
        output.flush();
        std::cout.rdbuf(previous_cout);
        std::clog.rdbuf(previous_clog);
    }
};

struct RandomStreams {
    explicit RandomStreams(int seed) {
        init_random(seed);
        rstream_vec.resize(1);
        init_random(seed + 1, false, &rstream_vec[0]);
        generator_vec.resize(1);
        generator_vec[0].seed(static_cast<unsigned int>(seed) + 2U);
    }

    ~RandomStreams() {
        finish_multi_rstreams();
        generator_vec.clear();
        finish_random();
    }
};

void simulate(const char* tree_file, const char* output_prefix, const char* model,
              const char* sequence_type, int length, int seed,
              const char* indel, const char* indel_size, bool no_unaligned,
              pf_alisim_output callback, void* context) {
    // Params borrows these strings; retain storage until the next call.
    // The Rust wrapper serializes all calls.
    static std::vector<std::string> arguments;
    arguments = {
        "pf-sim-alisim", "--alisim", output_prefix, "-t", tree_file,
        "-m", model, "--seqtype", sequence_type, "--length", std::to_string(length),
        "-af", "fasta", "--threads", "1", "--seed", std::to_string(seed)
    };
    if (indel) {
        arguments.emplace_back("--indel");
        arguments.emplace_back(indel);
    }
    if (indel_size) {
        arguments.emplace_back("--indel-size");
        arguments.emplace_back(indel_size);
    }
    if (no_unaligned) arguments.emplace_back("--no-unaligned");

    std::vector<char*> argv;
    for (auto& argument : arguments) argv.push_back(argument.data());
    const auto argc = static_cast<int>(argv.size());
    argv.push_back(nullptr);

    auto& params = Params::getInstance();
    parseArg(argc, argv.data(), params);
    params.SSE = LK_386;
    params.suppress_output_flags |= OUT_LOG;
    static const bool initialized = [] {
        read_distributions();
        return true;
    }();
    (void)initialized;

    RandomStreams random_streams(seed);
    pf::MemorySession memory;
    Checkpoint checkpoint;
    runAliSim(params, &checkpoint);
    if (pf::raw_enabled()) {
        for (const auto& row : pf::raw_state().rows) {
            if (callback(row.states.data(), row.states.size(), 0, context) != 0 ||
                callback(reinterpret_cast<const uint8_t*>(row.name.data()), row.name.size(), 1, context) != 0)
                throw std::runtime_error("Rust could not accept a raw row");
        }
        const uint8_t export_unaligned = params.alisim_insertion_ratio + params.alisim_deletion_ratio > 0
            && !params.alisim_no_export_sequence_wo_gaps;
        if (callback(&export_unaligned, 1, 3, context) != 0)
            throw std::runtime_error("Rust could not accept numeric output metadata");
        const auto& alphabet = pf::raw_state().alphabet;
        if (callback(alphabet.data(), alphabet.size(), 2, context) != 0)
            throw std::runtime_error("Rust could not accept the raw alphabet");
        return;
    }
    const auto fasta = std::string(output_prefix) + ".fa";
    const auto aligned = pf::memory_files().find(fasta);
    if (aligned == pf::memory_files().end() || aligned->second.empty()) {
        throw std::runtime_error("AliSim did not produce an in-memory FASTA alignment");
    }
    if (callback(reinterpret_cast<const uint8_t*>(aligned->second.data()),
                 aligned->second.size(), 0, context) != 0) {
        throw std::runtime_error("Rust could not accept the AliSim alignment");
    }
    const auto unaligned = pf::memory_files().find(std::string(output_prefix) + ".unaligned.fa");
    if (unaligned != pf::memory_files().end() &&
        callback(reinterpret_cast<const uint8_t*>(unaligned->second.data()),
                 unaligned->second.size(), 1, context) != 0) {
        throw std::runtime_error("Rust could not accept the AliSim gap-free sequences");
    }
}

}  // namespace

extern "C" int pf_alisim_simulate(const char* tree_file, const pf_tree_edge* edges,
                                 size_t edge_count, uint32_t tips, const char* output_prefix,
                                 const char* model, const char* sequence_type,
                                 int32_t length, int32_t seed,
                                 const char* indel, const char* indel_size, int no_unaligned,
                                 const char* log_path, int output_mode,
                                 pf_alisim_output callback, void* context,
                                 char* error, size_t error_capacity) {
    if (!error || error_capacity == 0) return 1;
    error[0] = '\0';
    auto fail = [&](const char* message) {
        const auto size = std::min(std::strlen(message), error_capacity - 1);
        std::memcpy(error, message, size);
        error[size] = '\0';
        return 1;
    };

    try {
        if (!tree_file || !output_prefix || !model || !sequence_type || !callback ||
            length <= 0 || seed <= 0 || seed >= INT32_MAX) {
            throw std::invalid_argument("Invalid AliSim bridge arguments");
        }
        OutputRedirect output(log_path);
        if (output_mode != 0 && output_mode != 1) throw std::invalid_argument("Invalid output mode");
        pf::RawSession raw(output_mode == 1);
        pf::TreeSession input_tree(edges, edge_count, tips);
        simulate(tree_file, output_prefix, model, sequence_type, length, seed,
                 indel, indel_size, no_unaligned != 0, callback, context);
        output.finish();
        return 0;
    } catch (const std::exception& exception) {
        return fail(exception.what());
    } catch (...) {
        return fail("Unknown C++ exception from AliSim");
    }
}
