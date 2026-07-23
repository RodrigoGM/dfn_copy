#include "scratch_matrix.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <stdexcept>
#include <cstring>

ScratchMatrix::ScratchMatrix(int fd, std::string path, size_t num_cells, size_t num_bins)
    : fd_(fd), path_(std::move(path)), num_cells_(num_cells), num_bins_(num_bins) {}

ScratchMatrix::ScratchMatrix(ScratchMatrix&& other) noexcept
    : fd_(other.fd_), path_(std::move(other.path_)),
      num_cells_(other.num_cells_), num_bins_(other.num_bins_) {
    other.fd_ = -1;
}

ScratchMatrix::~ScratchMatrix() {
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
}

ScratchMatrix ScratchMatrix::create(const std::string& path, size_t num_cells, size_t num_bins) {
    int fd = ::open(path.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0644);
    if (fd < 0) {
        throw std::runtime_error("cannot create scratch file: " + path +
                                  " (" + std::strerror(errno) + ")");
    }
    off_t total_bytes = static_cast<off_t>(num_cells) * static_cast<off_t>(num_bins) *
                         static_cast<off_t>(sizeof(float));
    if (::ftruncate(fd, total_bytes) != 0) {
        std::string msg = "cannot size scratch file: " + path + " (" + std::strerror(errno) + ")";
        ::close(fd);
        throw std::runtime_error(msg);
    }
    return ScratchMatrix(fd, path, num_cells, num_bins);
}

void ScratchMatrix::write_cell(size_t cell_index, const std::vector<double>& values) {
    if (values.size() != num_bins_) {
        throw std::runtime_error("scratch matrix: write_cell got " +
                                  std::to_string(values.size()) + " values, expected " +
                                  std::to_string(num_bins_));
    }
    if (cell_index >= num_cells_) {
        throw std::runtime_error("scratch matrix: cell index " + std::to_string(cell_index) +
                                  " out of bounds (num_cells=" + std::to_string(num_cells_) + ")");
    }
    std::vector<float> buf(num_bins_);
    for (size_t i = 0; i < num_bins_; ++i) buf[i] = static_cast<float>(values[i]);

    off_t offset = static_cast<off_t>(cell_index) * static_cast<off_t>(num_bins_) *
                   static_cast<off_t>(sizeof(float));
    size_t bytes_to_write = num_bins_ * sizeof(float);
    ssize_t written = ::pwrite(fd_, buf.data(), bytes_to_write, offset);
    if (written < 0 || static_cast<size_t>(written) != bytes_to_write) {
        throw std::runtime_error("scratch matrix: short/failed write for cell " +
                                  std::to_string(cell_index) + " (" + path_ + ")");
    }
}

std::vector<double> ScratchMatrix::read_cell(size_t cell_index) const {
    if (cell_index >= num_cells_) {
        throw std::runtime_error("scratch matrix: cell index " + std::to_string(cell_index) +
                                  " out of bounds (num_cells=" + std::to_string(num_cells_) + ")");
    }
    std::vector<float> buf(num_bins_);
    off_t offset = static_cast<off_t>(cell_index) * static_cast<off_t>(num_bins_) *
                   static_cast<off_t>(sizeof(float));
    size_t bytes_to_read = num_bins_ * sizeof(float);
    ssize_t read_n = ::pread(fd_, buf.data(), bytes_to_read, offset);
    if (read_n < 0 || static_cast<size_t>(read_n) != bytes_to_read) {
        throw std::runtime_error("scratch matrix: short/failed read for cell " +
                                  std::to_string(cell_index) + " (" + path_ + ")");
    }
    std::vector<double> out(num_bins_);
    for (size_t i = 0; i < num_bins_; ++i) out[i] = static_cast<double>(buf[i]);
    return out;
}

std::vector<float> ScratchMatrix::read_all() const {
    size_t total = num_cells_ * num_bins_;
    std::vector<float> out(total);
    size_t total_bytes = total * sizeof(float);
    size_t done = 0;
    while (done < total_bytes) {
        ssize_t n = ::pread(fd_, reinterpret_cast<char*>(out.data()) + done,
                            total_bytes - done, static_cast<off_t>(done));
        if (n <= 0) {
            throw std::runtime_error("scratch matrix: short/failed read_all (" + path_ + ")");
        }
        done += static_cast<size_t>(n);
    }
    return out;
}

void ScratchMatrix::close() {
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
}
