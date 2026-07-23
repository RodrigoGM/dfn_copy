#pragma once
#include <string>
#include <vector>
#include <mutex>
#include <fstream>
#include "segment_cell.hpp"

// Writes an IGV/DNAcopy-format .seg file (plain text, not gzipped --
// row-per-segment, not a matrix). Thread-safe: append_rows() is called
// once per completed cell from worker threads, guarded by an internal
// mutex.
class SegFileWriter {
public:
    explicit SegFileWriter(const std::string& path); // throws on open failure
    ~SegFileWriter();

    SegFileWriter(const SegFileWriter&) = delete;
    SegFileWriter& operator=(const SegFileWriter&) = delete;

    void append_rows(const std::vector<SegRow>& rows);

    // Idempotent. Safe to call more than once.
    void close();

private:
    std::ofstream out_;
    std::mutex mutex_;
};
