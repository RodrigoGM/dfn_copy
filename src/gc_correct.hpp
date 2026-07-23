#pragma once
#include <vector>

struct GcCorrectedCell {
    std::vector<double> gc_corrected;  // raw-count-scale, rescaled to this cell's median
    std::vector<double> lowess_ratio;  // dimensionless, mean approx 1
};

// GC-corrects one cell's bin-count vector against the bins' GC values.
// counts.size() must equal gc.size() -- throws std::runtime_error
// otherwise. If median(counts) == 0, both outputs equal counts unchanged
// (matching correct_gc.py's existing safeguard).
GcCorrectedCell gc_correct_cell(const std::vector<double>& counts,
                                 const std::vector<double>& gc);
