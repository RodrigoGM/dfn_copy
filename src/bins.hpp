#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <optional>
#include "cli_args.hpp"

struct Bin {
    std::string chrom;
    int64_t start;
    int64_t end;
    double gc;
};

class BinIndex {
public:
    static BinIndex load(const std::string& path);

    size_t size() const { return bins_.size(); }
    std::string name(size_t i) const;
    std::optional<size_t> find(const std::string& chrom, int64_t coord) const;
    const std::vector<Bin>& bins() const { return bins_; }

private:
    std::vector<Bin> bins_;
    std::unordered_map<std::string, std::vector<std::pair<int64_t, size_t>>> by_chrom_;
};

int64_t compute_bin_coordinate(int64_t align_start, int64_t align_end, PositionMode mode);

std::string validate_chrom_naming(const std::vector<std::string>& bam_contigs,
                                   const std::vector<std::string>& bin_chroms);
