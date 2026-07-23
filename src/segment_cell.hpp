#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "bin_gc.hpp"
#include "cbs_args.hpp"

struct SegRow {
    std::string barcode;
    std::string chrom;
    int64_t start;
    int64_t end;
    int num_mark;
    double seg_mean;
};

struct CellSegmentation {
    std::vector<double> segmented_lowess_ratio; // same length/order as bins
    std::vector<SegRow> seg_rows;
};

// Groups lowess_ratio by chromosome (first-seen order, per bins' chrom
// column), segments each chromosome independently via the vendored CBS
// core with centering forced off, and expands the result back to a
// bins-length vector plus this cell's .seg rows. bins.size() must equal
// lowess_ratio.size().
CellSegmentation segment_cell(const std::vector<BinGc>& bins,
                               const std::vector<double>& lowess_ratio,
                               const std::string& barcode,
                               const CbsArgs& cbs_args);
