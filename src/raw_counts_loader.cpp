#include "raw_counts_loader.hpp"
#include <zlib.h>
#include <sstream>
#include <stdexcept>
#include <limits>

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

// Reads one '\n'-terminated line from a gzFile. Returns false at EOF.
bool gz_getline(gzFile f, std::string& out) {
    out.clear();
    char buf[4096];
    bool got_any = false;
    while (gzgets(f, buf, sizeof(buf)) != nullptr) {
        got_any = true;
        out += buf;
        if (!out.empty() && out.back() == '\n') {
            out.pop_back();
            return true;
        }
        if (gzeof(f)) break;
    }
    return got_any;
}

int32_t parse_i32(const std::string& s, const std::string& context) {
    try {
        size_t consumed = 0;
        long v = std::stol(s, &consumed);
        if (consumed != s.size()) throw std::invalid_argument("trailing characters");
        if (v < static_cast<long>(std::numeric_limits<int32_t>::min()) ||
            v > static_cast<long>(std::numeric_limits<int32_t>::max())) {
            throw std::runtime_error("raw counts file: value out of int32 range for " +
                                      context + ": '" + s + "'");
        }
        return static_cast<int32_t>(v);
    } catch (const std::runtime_error&) {
        throw;
    } catch (const std::exception&) {
        throw std::runtime_error("raw counts file: expected an integer for " +
                                  context + ", got '" + s + "'");
    }
}

} // namespace

RawCountsMatrix load_raw_counts(const std::string& path,
                                 size_t expected_num_bins,
                                 const std::vector<std::string>& expected_bin_names) {
    gzFile in = gzopen(path.c_str(), "rb");
    if (!in) {
        throw std::runtime_error("cannot open raw counts file: " + path);
    }

    std::string header_line;
    if (!gz_getline(in, header_line)) {
        int errnum = 0;
        const char* err_msg = gzerror(in, &errnum);
        if (errnum != Z_OK && errnum != Z_STREAM_END) {
            std::string msg = "raw counts file is corrupted or not a valid gzip file: " +
                               path + " (" + err_msg + ")";
            gzclose(in);
            throw std::runtime_error(msg);
        }
        gzclose(in);
        throw std::runtime_error("raw counts file is empty: " + path);
    }
    std::vector<std::string> header = split_tab(header_line);
    if (header.empty() || header[0] != "bin") {
        gzclose(in);
        throw std::runtime_error("raw counts file: expected header starting with "
                                  "'bin' (" + path + ")");
    }

    RawCountsMatrix m;
    m.num_cells = header.size() - 1;
    m.barcodes.assign(header.begin() + 1, header.end());
    m.num_bins = expected_num_bins;

    if (expected_bin_names.size() != expected_num_bins) {
        gzclose(in);
        throw std::runtime_error("raw counts loader: expected_bin_names size (" +
            std::to_string(expected_bin_names.size()) + ") does not match expected_num_bins (" +
            std::to_string(expected_num_bins) + ") -- this is a caller contract violation");
    }

    m.data.assign(m.num_bins * m.num_cells, 0);

    std::string line;
    size_t bin_row = 0;
    size_t line_no = 1;
    while (gz_getline(in, line)) {
        ++line_no;
        if (line.empty()) continue;
        if (bin_row >= expected_num_bins) {
            gzclose(in);
            throw std::runtime_error("raw counts file: has more bin rows than "
                                      "the bin file (" + path + ")");
        }
        std::vector<std::string> fields = split_tab(line);
        if (fields.size() != m.num_cells + 1) {
            gzclose(in);
            throw std::runtime_error("raw counts file: line " + std::to_string(line_no) +
                                      " has " + std::to_string(fields.size() - 1) +
                                      " value column(s), expected " +
                                      std::to_string(m.num_cells) + " (" + path + ")");
        }
        if (fields[0] != expected_bin_names[bin_row]) {
            gzclose(in);
            throw std::runtime_error("raw counts file: bin row " + std::to_string(bin_row) +
                                      " is '" + fields[0] + "', expected '" +
                                      expected_bin_names[bin_row] +
                                      "' -- --counts and --bins must be the matching pair "
                                      "from the same dfn_copy run (" + path + ")");
        }
        for (size_t c = 0; c < m.num_cells; ++c) {
            try {
                m.data[bin_row * m.num_cells + c] =
                    parse_i32(fields[c + 1], "value at line " + std::to_string(line_no));
            } catch (const std::runtime_error&) {
                gzclose(in);
                throw;
            }
        }
        ++bin_row;
    }
    gzclose(in);

    if (bin_row != expected_num_bins) {
        throw std::runtime_error("raw counts file: has " + std::to_string(bin_row) +
                                  " bin row(s), expected " + std::to_string(expected_num_bins) +
                                  " to match --bins (" + path + ")");
    }

    return m;
}
