#include "lowess.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>

namespace {

// One tri-cube-weighted local linear fit at every point, given the
// current per-point robustness weights (all 1.0 on the initial pass).
// order/xs are x's ascending-sorted permutation/values, precomputed once
// by the caller and reused across passes (x doesn't change between
// passes -- only the robustness weights do).
std::vector<double> weighted_fit_pass(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<size_t>& order,
                                       const std::vector<double>& xs,
                                       const std::vector<double>& robustness,
                                       size_t k) {
    size_t n = x.size();
    std::vector<double> fitted(n, 0.0);

    // inv_order[order[t]] = t : original index -> position in sorted order.
    std::vector<size_t> inv_order(n);
    for (size_t t = 0; t < n; ++t) inv_order[order[t]] = t;

    for (size_t t = 0; t < n; ++t) {
        size_t orig_i = order[t];
        double xi = xs[t];

        // Find the k nearest neighbors (by x-distance) around sorted
        // position t, via two-pointer expansion.
        size_t lo = t, hi = t;
        while (hi - lo + 1 < k) {
            bool can_left = lo > 0;
            bool can_right = hi + 1 < n;
            if (can_left && can_right) {
                double dl = xi - xs[lo - 1];
                double dr = xs[hi + 1] - xi;
                if (dl <= dr) --lo; else ++hi;
            } else if (can_left) {
                --lo;
            } else if (can_right) {
                ++hi;
            } else {
                break; // window already covers all n points
            }
        }

        double d_max = std::max(xi - xs[lo], xs[hi] - xi);

        double S0 = 0.0, S1 = 0.0, S2 = 0.0, Sy = 0.0, Sxy = 0.0;
        for (size_t p = lo; p <= hi; ++p) {
            size_t orig_j = order[p];
            double xj = xs[p];
            double dist = std::fabs(xj - xi);
            double tricube = 1.0;
            if (d_max > 0.0) {
                double u = dist / d_max;
                if (u >= 1.0) tricube = 0.0;
                else {
                    double v = 1.0 - u * u * u;
                    tricube = v * v * v;
                }
            }
            double w = tricube * robustness[orig_j];
            double yj = y[orig_j];
            S0 += w;
            S1 += w * xj;
            S2 += w * xj * xj;
            Sy += w * yj;
            Sxy += w * xj * yj;
        }

        double denom = S0 * S2 - S1 * S1;
        double a, b;
        if (S0 <= 0.0) {
            // Degenerate: no positive weight anywhere in the window
            // (shouldn't happen with a well-formed robustness vector, but
            // guarded rather than dividing by zero). Fall back to the raw
            // value itself so the point simply isn't smoothed this pass.
            fitted[orig_i] = y[orig_i];
            continue;
        }
        if (std::fabs(denom) < 1e-12 * std::max(1.0, S0 * S0)) {
            // Colinear/degenerate window (e.g. every x in range is
            // identical) -- fall back to the weighted mean.
            a = Sy / S0;
            b = 0.0;
        } else {
            b = (S0 * Sxy - S1 * Sy) / denom;
            a = (Sy - b * S1) / S0;
        }
        fitted[orig_i] = a + b * xi;
        (void)inv_order; // computed for clarity/future use; not otherwise read
    }

    return fitted;
}

double median_abs(const std::vector<double>& v) {
    std::vector<double> a(v.size());
    for (size_t i = 0; i < v.size(); ++i) a[i] = std::fabs(v[i]);
    std::sort(a.begin(), a.end());
    size_t n = a.size();
    if (n == 0) return 0.0;
    if (n % 2 == 1) return a[n / 2];
    return 0.5 * (a[n / 2 - 1] + a[n / 2]);
}

} // namespace

std::vector<double> lowess(const std::vector<double>& y,
                            const std::vector<double>& x,
                            double frac,
                            int iterations) {
    if (x.size() != y.size()) {
        throw std::runtime_error("lowess: x and y must be the same size");
    }
    size_t n = x.size();
    if (n == 0) {
        throw std::runtime_error("lowess: x/y must be non-empty");
    }

    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return x[a] < x[b]; });

    std::vector<double> xs(n);
    for (size_t t = 0; t < n; ++t) xs[t] = x[order[t]];

    size_t k = static_cast<size_t>(std::ceil(frac * static_cast<double>(n)));
    if (k < 2) k = 2;
    if (k > n) k = n;

    std::vector<double> robustness(n, 1.0);
    std::vector<double> fitted = weighted_fit_pass(x, y, order, xs, robustness, k);

    for (int pass = 0; pass < iterations; ++pass) {
        std::vector<double> residual(n);
        for (size_t i = 0; i < n; ++i) residual[i] = y[i] - fitted[i];
        double s = median_abs(residual);
        if (s <= 0.0) {
            // Fit is already exact everywhere; further reweighting can't
            // change anything (bisquare of all-zero residuals is
            // undefined via division by s -- but there's nothing left to
            // robustify against).
            break;
        }
        for (size_t i = 0; i < n; ++i) {
            double u = residual[i] / (6.0 * s);
            if (u * u >= 1.0) {
                robustness[i] = 0.0;
            } else {
                double v = 1.0 - u * u;
                robustness[i] = v * v;
            }
        }
        fitted = weighted_fit_pass(x, y, order, xs, robustness, k);
    }

    return fitted;
}
