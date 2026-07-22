#pragma once
#include <vector>

// Locally weighted scatterplot smoothing (Cleveland 1979), matching
// statsmodels.nonparametric.smoothers_lowess.lowess's default behavior:
// frac fraction of points used as the local window (tri-cube-weighted),
// `iterations` residual-based bisquare reweighting passes performed after
// an initial fit, local *linear* fit at each point. Returns fitted values
// in the same order as the input x/y (not sorted, not resampled onto a
// grid). x and y must be the same size and non-empty.
std::vector<double> lowess(const std::vector<double>& y,
                            const std::vector<double>& x,
                            double frac = 2.0 / 3.0,
                            int iterations = 3);
