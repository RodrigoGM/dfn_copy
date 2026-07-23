#pragma once
#include <vector>
#include <tuple>

struct MaxStat { double stat; int k; }; // 1-CP split after k (s..k | k+1..e)
struct MaxStat2 { double stat; int k1; int k2; }; // 2-CP (s..k1 | k1+1..k2 | k2+1..e)

MaxStat max_cusum_1cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e);
MaxStat2 max_cusum_2cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e, int min_seg_len);

inline double seg_mean(const std::vector<double>& csum, int s, int e){ return (csum[e+1]-csum[s]) / (double)(e-s+1); }
void build_csum(const std::vector<double>& x, std::vector<double>& csum);
