#include "bin_gc.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

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

std::vector<BinGc> load_bin_gc(const std::string& path) {
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

    std::vector<BinGc> bins;
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
        BinGc b;
        b.chrom = fields[i_chrom];
        std::string ctx = "line " + std::to_string(line_no);
        b.start = parse_i64(fields[i_start], "start at " + ctx);
        b.end = parse_i64(fields[i_end], "end at " + ctx);
        b.gc = parse_f64(fields[i_gc], "gc at " + ctx);
        bins.push_back(b);
    }

    return bins;
}

std::string bin_name(const BinGc& b) {
    return b.chrom + ":" + std::to_string(b.start) + ":" + std::to_string(b.end);
}
