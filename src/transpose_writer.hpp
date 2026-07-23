#pragma once
#include <string>
#include <vector>

// Writes a cell-major flat float32 array (as loaded from
// ScratchMatrix::read_all() -- row-major by cell,
// scratch_data[cell * num_bins + bin]) out as a bins x cells gzip TSV,
// the same convention dfn_copy's raw_counts.txt.gz uses. Throws
// std::runtime_error on any size mismatch or an unwritable path.
void write_transposed_gz(const std::string& path,
                          const std::vector<float>& scratch_data,
                          size_t num_cells,
                          size_t num_bins,
                          const std::vector<std::string>& bin_names,
                          const std::vector<std::string>& barcode_names);
