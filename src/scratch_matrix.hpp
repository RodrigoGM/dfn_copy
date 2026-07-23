#pragma once
#include <cstddef>
#include <string>
#include <vector>

// A flat binary scratch file for one cell-major output matrix: fixed
// record size (num_bins float32s per cell), opened once, written to
// concurrently by many threads via positioned I/O (pwrite/pread) --
// no locking needed since each cell's byte range is disjoint from every
// other cell's.
class ScratchMatrix {
public:
    static ScratchMatrix create(const std::string& path, size_t num_cells, size_t num_bins);

    ScratchMatrix(const ScratchMatrix&) = delete;
    ScratchMatrix& operator=(const ScratchMatrix&) = delete;
    ScratchMatrix(ScratchMatrix&& other) noexcept;
    ScratchMatrix& operator=(ScratchMatrix&&) = delete;
    ~ScratchMatrix();

    // values.size() must equal num_bins(). Thread-safe across distinct
    // cell_index values; not safe to call twice concurrently for the
    // *same* cell_index (not needed -- each cell is written exactly once).
    void write_cell(size_t cell_index, const std::vector<double>& values);

    // Reads one cell back as double (upcast from the file's float32
    // storage) -- used for round-trip testing, not the hot path.
    std::vector<double> read_cell(size_t cell_index) const;

    // Reads the entire file into memory as a flat, row-major-by-cell
    // float32 vector (size num_cells() * num_bins()) -- what Task 12's
    // transpose-on-exit uses.
    std::vector<float> read_all() const;

    size_t num_cells() const { return num_cells_; }
    size_t num_bins() const { return num_bins_; }
    const std::string& path() const { return path_; }

    // Idempotent. Safe to call more than once.
    void close();

private:
    ScratchMatrix(int fd, std::string path, size_t num_cells, size_t num_bins);
    int fd_;
    std::string path_;
    size_t num_cells_;
    size_t num_bins_;
};
