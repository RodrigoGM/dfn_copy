#include "test_framework.hpp"
#include "lowess.hpp"
#include <cmath>

namespace {
double mean(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

double pearson(const std::vector<double>& a, const std::vector<double>& b) {
    double ma = mean(a), mb = mean(b);
    double num = 0.0, da = 0.0, db = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double xa = a[i] - ma, xb = b[i] - mb;
        num += xa * xb;
        da += xa * xa;
        db += xb * xb;
    }
    if (da <= 0.0 || db <= 0.0) return 0.0;
    return num / std::sqrt(da * db);
}
}

void test_recovers_a_clean_linear_trend() {
    // y = 100 + 200*x, x evenly spaced in [0,1]. A perfectly linear trend
    // should be recovered almost exactly by a local *linear* fit.
    const int n = 60;
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = static_cast<double>(i) / (n - 1);
        y[i] = 100.0 + 200.0 * x[i];
    }
    auto fitted = lowess(y, x);
    ASSERT_EQ(fitted.size(), static_cast<size_t>(n));
    for (int i = 5; i < n - 5; ++i) { // skip boundary points, more edge bias there
        ASSERT_NEAR(fitted[i], y[i], 5.0);
    }
}

void test_flattens_correlation_with_gc() {
    // Inject a strong linear GC trend into otherwise-flat counts, plus a
    // little per-point noise (deterministic pseudo-noise, not random --
    // keeps the test reproducible). Correlation(counts, gc) should be
    // high; correlation(counts - fitted, gc) should collapse toward 0,
    // mirroring correct_gc.py's own existing validation style.
    const int n = 80;
    std::vector<double> gc(n), counts(n);
    for (int i = 0; i < n; ++i) {
        gc[i] = static_cast<double>(i) / (n - 1);
        double noise = 3.0 * std::sin(static_cast<double>(i) * 0.7);
        counts[i] = 50.0 + 300.0 * gc[i] + noise;
    }
    double corr_before = pearson(counts, gc);
    ASSERT_TRUE(corr_before > 0.9);

    auto fitted = lowess(counts, gc);
    std::vector<double> residual(n);
    for (int i = 0; i < n; ++i) residual[i] = counts[i] - fitted[i];
    double corr_after = std::fabs(pearson(residual, gc));
    ASSERT_TRUE(corr_after < 0.3);
}

void test_handles_unsorted_x_and_returns_original_order() {
    // x is NOT pre-sorted -- lowess() must sort internally and return
    // fitted values back in the caller's original order.
    std::vector<double> x = {0.8, 0.1, 0.5, 0.3, 0.9, 0.2, 0.6, 0.4, 0.7, 0.0};
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) y[i] = 10.0 + 20.0 * x[i];

    auto fitted = lowess(y, x);
    ASSERT_EQ(fitted.size(), x.size());
    // x[9] = 0.0 is the smallest value; its fit should be near y[9] = 10.0.
    ASSERT_NEAR(fitted[9], 10.0, 5.0);
    // x[4] = 0.9 is the largest value; its fit should be near y[4] = 28.0.
    ASSERT_NEAR(fitted[4], 28.0, 5.0);
}

int main() {
    test_recovers_a_clean_linear_trend();
    test_flattens_correlation_with_gc();
    test_handles_unsorted_x_and_returns_original_order();
    TEST_REPORT();
}
