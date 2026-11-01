#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Calls must be serialized. Strings are borrowed and NUL-terminated.
// error addresses error_capacity writable bytes (capacity > 0).
// Returns 0 on success, 1 on a caught error; error is always NUL-terminated.
// Some errors inside IQ-TREE still terminate the process.
// log_path is NULL to discard informational output, or a borrowed path to
// append it to. Native errors remain on stderr.
// indel and indel_size are optional borrowed strings passed verbatim to AliSim.
// callback copies the aligned (kind=0) and optional gap-free (kind=1) FASTA
// bytes before returning. Pointers are borrowed only for the callback duration.
// A nonzero callback result aborts with an error. No MSA files are written.
// output_mode=0 is reference FASTA. Mode 1 is numeric output: row bytes, names,
// alphabet and a one-byte export-unaligned flag use callback kinds 0,1,2,3.
// Node IDs 0..tips are T1..Tn; other IDs are internal.
// Edge count 2*tips-3 means unrooted; 2*tips-2 means rooted (binary).
// The Newick display origin is node tips. Edges retain Rust's neighbor order.
typedef struct pf_tree_edge { uint32_t a; uint32_t b; double length; } pf_tree_edge;
// edges=NULL selects tree_file input. Otherwise tree_file is only a log label;
// edges are borrowed for the call and no tree file is opened.
typedef int (*pf_alisim_output)(const uint8_t* data, size_t size, int kind, void* context);
int pf_alisim_simulate(const char* tree_file, const pf_tree_edge* edges,
                      size_t edge_count, uint32_t tips, const char* output_prefix,
                      const char* model, const char* sequence_type,
                      int32_t length, int32_t seed,
                      const char* indel, const char* indel_size, int no_unaligned,
                      const char* log_path, int output_mode,
                      pf_alisim_output callback, void* context,
                      char* error, size_t error_capacity);

#ifdef __cplusplus
}
#endif
