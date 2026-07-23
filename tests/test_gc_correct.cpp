#include "test_framework.hpp"
#include "gc_correct.hpp"
#include <cmath>

void test_ratio_has_mean_near_one_and_corrected_median_matches_raw() {
    const int n = 60;
    std::vector<double> gc(n), counts(n);
    for (int i = 0; i < n; ++i) {
        gc[i] = static_cast<double>(i) / (n - 1);
        counts[i] = 80.0 + 100.0 * gc[i]; // clean GC trend, no noise
    }

    auto result = gc_correct_cell(counts, gc);
    ASSERT_EQ(result.gc_corrected.size(), static_cast<size_t>(n));
    ASSERT_EQ(result.lowess_ratio.size(), static_cast<size_t>(n));

    double ratio_mean = 0.0;
    for (double r : result.lowess_ratio) ratio_mean += r;
    ratio_mean /= n;
    ASSERT_NEAR(ratio_mean, 1.0, 0.1);

    std::vector<double> raw_sorted = counts, corrected_sorted = result.gc_corrected;
    std::sort(raw_sorted.begin(), raw_sorted.end());
    std::sort(corrected_sorted.begin(), corrected_sorted.end());
    double raw_median = raw_sorted[n / 2];
    double corrected_median = corrected_sorted[n / 2];
    ASSERT_NEAR(corrected_median, raw_median, raw_median * 0.15);
}

void test_zero_median_falls_back_to_raw_counts() {
    std::vector<double> gc = {0.1, 0.2, 0.3, 0.4, 0.5};
    std::vector<double> counts = {0.0, 0.0, 0.0, 0.0, 0.0}; // median is 0
    auto result = gc_correct_cell(counts, gc);
    for (size_t i = 0; i < counts.size(); ++i) {
        ASSERT_NEAR(result.gc_corrected[i], counts[i], 1e-9);
        ASSERT_NEAR(result.lowess_ratio[i], counts[i], 1e-9);
    }
}

void test_mismatched_sizes_throws() {
    std::vector<double> gc = {0.1, 0.2};
    std::vector<double> counts = {1.0, 2.0, 3.0};
    bool threw = false;
    try {
        gc_correct_cell(counts, gc);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

void test_empty_input_throws() {
    std::vector<double> gc;
    std::vector<double> counts;
    bool threw = false;
    try {
        gc_correct_cell(counts, gc);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

int main() {
    test_ratio_has_mean_near_one_and_corrected_median_matches_raw();
    test_zero_median_falls_back_to_raw_counts();
    test_mismatched_sizes_throws();
    test_empty_input_throws();
    TEST_REPORT();
}
