#include "gc_correct.hpp"
#include "lowess.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

GcCorrectedCell gc_correct_cell(const std::vector<double>& counts,
                                 const std::vector<double>& gc) {
    if (counts.size() != gc.size()) {
        throw std::runtime_error("gc_correct_cell: counts and gc must be the same size");
    }
    size_t n = counts.size();
    GcCorrectedCell result;
    result.gc_corrected.resize(n);
    result.lowess_ratio.resize(n);

    std::vector<double> sorted_counts = counts;
    std::sort(sorted_counts.begin(), sorted_counts.end());
    double median = (n % 2 == 1)
        ? sorted_counts[n / 2]
        : 0.5 * (sorted_counts[n / 2 - 1] + sorted_counts[n / 2]);

    if (median == 0.0) {
        result.gc_corrected = counts;
        result.lowess_ratio = counts;
        return result;
    }

    std::vector<double> fitted = lowess(counts, gc);
    const double eps = std::numeric_limits<double>::epsilon();
    for (size_t i = 0; i < n; ++i) {
        double f = fitted[i];
        if (f <= 0.0) f = eps;
        result.lowess_ratio[i] = counts[i] / f;
        result.gc_corrected[i] = counts[i] * (median / f);
    }

    return result;
}
