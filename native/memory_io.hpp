#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

// The build adapts only AliSim's uncompressed sequence output streams and its
// indel temporary-data reader. The Rust wrapper serializes native calls.
namespace pf {

inline std::unordered_map<std::string, std::string>& memory_files() {
    static std::unordered_map<std::string, std::string> files;
    return files;
}

struct MemorySession {
    MemorySession() { memory_files().clear(); }
    ~MemorySession() { memory_files().clear(); }
};

class MemoryOutputFile : public std::ostringstream {
    std::string path;
    bool closed = false;

public:
    explicit MemoryOutputFile(const char* filename,
                              std::ios_base::openmode mode = std::ios_base::out)
        : std::ostringstream(mode | std::ios_base::out), path(filename) {
        if (mode & std::ios_base::app) {
            str(memory_files()[path]);
            seekp(0, std::ios_base::end);
        }
    }

    void close() {
        if (!closed) {
            memory_files()[path] = str();
            closed = true;
        }
    }
};

class MemoryInputFile : public std::istringstream {
public:
    void open(const char* filename) {
        const auto found = memory_files().find(filename);
        if (found == memory_files().end()) {
            throw std::runtime_error("Missing AliSim in-memory temporary data");
        }
        str(found->second);
        clear();
    }

    void close() { str(""); }
};

}  // namespace pf
