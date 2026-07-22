#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct BinGc {
    std::string chrom;
    int64_t start;
    int64_t end;
    double gc;
};

// Loads a bins.tsv (header-driven columns: chrom, start, end, gc, any
// order, extra columns ignored) in file order. Throws std::runtime_error
// with a specific message on a missing column or a malformed row.
std::vector<BinGc> load_bin_gc(const std::string& path);

std::string bin_name(const BinGc& b);
