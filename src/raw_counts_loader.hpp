#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct RawCountsMatrix {
    size_t num_bins = 0;
    size_t num_cells = 0;
    std::vector<std::string> barcodes;
    std::vector<int32_t> data; // row-major by bin: data[bin * num_cells + cell]
};

// Loads <prefix>.raw_counts.txt.gz fully into memory. Validates the file's
// bin-row count and names against expected_num_bins/expected_bin_names
// (in order) -- throws std::runtime_error with a specific message on any
// mismatch, malformed row, or unreadable file, never returning a partial
// or silently-wrong matrix.
RawCountsMatrix load_raw_counts(const std::string& path,
                                 size_t expected_num_bins,
                                 const std::vector<std::string>& expected_bin_names);
