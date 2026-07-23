#include "segment_cell.hpp"
#include "cbs.hpp"
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

CellSegmentation segment_cell(const std::vector<BinGc>& bins,
                               const std::vector<double>& lowess_ratio,
                               const std::string& barcode,
                               const CbsArgs& cbs_args) {
    if (bins.size() != lowess_ratio.size()) {
        throw std::runtime_error("segment_cell: bins and lowess_ratio must be the same size");
    }

    CellSegmentation result;
    result.segmented_lowess_ratio.assign(bins.size(), 0.0);

    // Group original bin indices by chromosome, preserving first-seen
    // chromosome order (robust to a non-contiguous chrom block, though
    // bins.tsv is expected to already be chromosome-grouped in practice).
    std::vector<std::string> chrom_order;
    std::unordered_map<std::string, std::vector<size_t>> by_chrom;
    for (size_t i = 0; i < bins.size(); ++i) {
        const std::string& c = bins[i].chrom;
        if (by_chrom.find(c) == by_chrom.end()) chrom_order.push_back(c);
        by_chrom[c].push_back(i);
    }

    Args a;
    a.alpha = cbs_args.alpha;
    a.perms = cbs_args.perms;
    a.min_seg_len = cbs_args.min_seg_len;
    a.max_depth = cbs_args.max_depth;
    a.method = (cbs_args.cbs_method == CbsMethod::OneCp) ? "1cp" : "2cp";
    a.center = false; // forced off -- see design spec's "Related work: cbs+"

    std::hash<std::string> hasher;

    for (const auto& chrom : chrom_order) {
        const auto& indices = by_chrom.at(chrom);

        Series S;
        S.chrom = chrom;
        S.pos.reserve(indices.size());
        S.x.reserve(indices.size());
        for (size_t idx : indices) {
            S.pos.push_back(static_cast<long>(bins[idx].start));
            S.x.push_back(lowess_ratio[idx]);
        }

        // std::hash<std::string> is implementation-defined (its values may
        // differ across compilers/standard-library versions/platforms), so
        // this seed -- and therefore the exact permutation draws below -- is
        // reproducible only within a single build/toolchain, not bit-for-bit
        // across machines. Segmentation is a statistical permutation test, so
        // cross-platform results are statistically equivalent, not identical.
        uint64_t chrom_seed = cbs_args.seed ^ hasher(barcode) ^ hasher(chrom);
        std::mt19937_64 rng(chrom_seed);

        std::vector<Segment> segs;
        segment_chromosome(S, a, segs, rng);
        std::sort(segs.begin(), segs.end(),
                  [](const Segment& x, const Segment& y) { return x.s < y.s; });

        for (const auto& seg : segs) {
            for (int local = seg.s; local <= seg.e; ++local) {
                result.segmented_lowess_ratio[indices[local]] = seg.mean;
            }
            SegRow row;
            row.barcode = barcode;
            row.chrom = chrom;
            row.start = bins[indices[seg.s]].start;
            row.end = bins[indices[seg.e]].end;
            row.num_mark = seg.e - seg.s + 1;
            row.seg_mean = seg.mean;
            result.seg_rows.push_back(row);
        }
    }

    return result;
}
