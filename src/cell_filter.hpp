#pragma once
#include <cstdint>
#include <vector>
#include "raw_counts_loader.hpp"

struct FilteredCells {
    std::vector<size_t> surviving_indices; // original column indices, ascending
    int64_t max_column_sum = 0;
};

// Sums each cell's (column's) raw counts across all bins and keeps only
// cells whose sum is >= min_reads. max_column_sum is always the largest
// observed column sum, regardless of how many cells survive -- used for
// error reporting when the filter is too strict.
FilteredCells filter_cells(const RawCountsMatrix& m, int64_t min_reads);
