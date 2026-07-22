#include "bins.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>

namespace {

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }
    return fields;
}

int64_t parse_i64(const std::string& s, const std::string& context) {
    try {
        size_t consumed = 0;
        int64_t v = std::stoll(s, &consumed);
        if (consumed != s.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        throw std::runtime_error("bin file: expected an integer for " + context +
                                  ", got '" + s + "'");
    }
}

double parse_f64(const std::string& s, const std::string& context) {
    try {
        size_t consumed = 0;
        double v = std::stod(s, &consumed);
        if (consumed != s.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        throw std::runtime_error("bin file: expected a number for " + context +
                                  ", got '" + s + "'");
    }
}

} // namespace

BinIndex BinIndex::load(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot open bin file: " + path);
    }

    std::string header_line;
    if (!std::getline(in, header_line)) {
        throw std::runtime_error("bin file is empty: " + path);
    }
    std::vector<std::string> header = split_tab(header_line);

    auto col_index = [&](const std::string& name) -> int {
        for (size_t i = 0; i < header.size(); ++i) {
            if (header[i] == name) return static_cast<int>(i);
        }
        return -1;
    };

    int i_chrom = col_index("chrom");
    int i_start = col_index("start");
    int i_end = col_index("end");
    int i_gc = col_index("gc");

    std::vector<std::string> missing;
    if (i_chrom < 0) missing.push_back("chrom");
    if (i_start < 0) missing.push_back("start");
    if (i_end < 0) missing.push_back("end");
    if (i_gc < 0) missing.push_back("gc");
    if (!missing.empty()) {
        std::string msg = "bin file missing required column(s):";
        for (const auto& m : missing) msg += " " + m;
        throw std::runtime_error(msg + " (" + path + ")");
    }

    BinIndex idx;
    std::string line;
    size_t line_no = 1;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty()) continue;
        std::vector<std::string> fields = split_tab(line);
        size_t max_needed = static_cast<size_t>(
            std::max({i_chrom, i_start, i_end, i_gc})) + 1;
        if (fields.size() < max_needed) {
            throw std::runtime_error("bin file: line " + std::to_string(line_no) +
                                      " has too few columns (" + path + ")");
        }
        Bin bin;
        bin.chrom = fields[i_chrom];
        std::string ctx = "line " + std::to_string(line_no);
        bin.start = parse_i64(fields[i_start], "start at " + ctx);
        bin.end = parse_i64(fields[i_end], "end at " + ctx);
        bin.gc = parse_f64(fields[i_gc], "gc at " + ctx);

        size_t original_index = idx.bins_.size();
        idx.bins_.push_back(bin);
        idx.by_chrom_[bin.chrom].emplace_back(bin.start, original_index);
    }

    for (auto& [chrom, vec] : idx.by_chrom_) {
        std::sort(vec.begin(), vec.end());
    }

    return idx;
}

std::string BinIndex::name(size_t i) const {
    const Bin& b = bins_[i];
    return b.chrom + ":" + std::to_string(b.start) + ":" + std::to_string(b.end);
}

std::optional<size_t> BinIndex::find(const std::string& chrom, int64_t coord) const {
    auto it = by_chrom_.find(chrom);
    if (it == by_chrom_.end()) return std::nullopt;

    const auto& vec = it->second;
    // Find the last entry with start <= coord.
    auto upper = std::upper_bound(vec.begin(), vec.end(),
        std::make_pair(coord, static_cast<size_t>(-1)));
    if (upper == vec.begin()) return std::nullopt;
    --upper;
    size_t bin_idx = upper->second;
    const Bin& b = bins_[bin_idx];
    if (coord >= b.start && coord < b.end) {
        return bin_idx;
    }
    return std::nullopt;
}

int64_t compute_bin_coordinate(int64_t align_start, int64_t align_end, PositionMode mode) {
    switch (mode) {
        case PositionMode::Start:
            return align_start;
        case PositionMode::End:
            return align_end > align_start ? align_end - 1 : align_start;
        case PositionMode::Midpoint:
            return align_start + (align_end - align_start) / 2;
    }
    return align_start;
}

std::string validate_chrom_naming(const std::vector<std::string>& bam_contigs,
                                   const std::vector<std::string>& bin_chroms) {
    std::unordered_set<std::string> bam_set(bam_contigs.begin(), bam_contigs.end());
    std::vector<std::string> missing;
    for (const auto& c : bin_chroms) {
        if (bam_set.find(c) == bam_set.end()) missing.push_back(c);
    }
    if (missing.empty()) return "";

    auto has_chr_prefix = [](const std::string& s) {
        return s.size() > 3 && s.compare(0, 3, "chr") == 0;
    };

    size_t bam_chr_count = 0;
    for (const auto& c : bam_contigs) if (has_chr_prefix(c)) bam_chr_count++;
    size_t missing_chr_count = 0;
    for (const auto& c : missing) if (has_chr_prefix(c)) missing_chr_count++;

    bool bam_mostly_chr = !bam_contigs.empty() && bam_chr_count * 2 > bam_contigs.size();
    bool missing_mostly_chr = !missing.empty() && missing_chr_count * 2 > missing.size();

    std::ostringstream oss;
    if (bam_mostly_chr && !missing_mostly_chr) {
        oss << "Chromosome-naming mismatch: BAM header uses 'chr'-prefixed contig "
               "names (e.g. 'chr1'), but " << missing.size()
            << " bin chromosome(s) do not (e.g. '" << missing.front()
            << "'). Fix the bin file or BAM to use a consistent naming convention.";
    } else if (!bam_mostly_chr && missing_mostly_chr) {
        oss << "Chromosome-naming mismatch: bin file uses 'chr'-prefixed chromosome "
               "names (e.g. '" << missing.front()
            << "'), but the BAM header's contig names do not. Fix the bin file or "
               "BAM to use a consistent naming convention.";
    } else {
        oss << "Bin file references " << missing.size()
            << " chromosome(s) not found in the BAM header, e.g. '"
            << missing.front() << "'.";
    }
    return oss.str();
}
