#pragma once
#include <vector>
#include <string>
#include <unordered_map>

// Simple correction utilities: rolling-median detrend and linear covariate regression

void rolling_median_detrend(std::vector<double>& x, int window);

// If covariate columns are present: x <- x - (a + b1*c1 + b2*c2 + ...)
// Returns (intercept, coefficients) in out_beta for reference; silent if no covariates
void linear_covariate_regress(std::vector<double>& x, const std::unordered_map<std::string,std::vector<double>>& covars,
                              std::vector<double>& out_beta);
