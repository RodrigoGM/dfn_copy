#include "seg_file_writer.hpp"
#include <stdexcept>

SegFileWriter::SegFileWriter(const std::string& path) : path_(path) {
    out_.open(path, std::ios::out | std::ios::trunc);
    if (!out_.is_open()) {
        throw std::runtime_error("cannot open .seg file for writing: " + path);
    }
    out_ << "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n";
}

SegFileWriter::~SegFileWriter() {
    close();
}

void SegFileWriter::append_rows(const std::vector<SegRow>& rows) {
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto& row : rows) {
        out_ << row.barcode << '\t' << row.chrom << '\t' << row.start << '\t'
             << row.end << '\t' << row.num_mark << '\t' << row.seg_mean << '\n';
    }
    if (out_.fail()) {
        throw std::runtime_error("seg file writer: write failed (disk full or I/O error) for " +
                                  path_);
    }
}

void SegFileWriter::close() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (out_.is_open()) {
        out_.close();
    }
}
