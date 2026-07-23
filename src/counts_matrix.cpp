#include "counts_matrix.hpp"
#include <zlib.h>
#include <sstream>
#include <stdexcept>

CountsMatrix::CountsMatrix(size_t num_bins, size_t initial_barcodes)
    : num_bins_(num_bins),
      num_barcodes_(initial_barcodes),
      stride_(initial_barcodes == 0 ? 1 : initial_barcodes),
      data_(num_bins_ * stride_, 0) {}

void CountsMatrix::grow_to(size_t n) {
    if (n <= num_barcodes_) return;

    if (n > stride_) {
        size_t new_stride = stride_ == 0 ? 1 : stride_;
        while (new_stride < n) new_stride *= 2;

        std::vector<int32_t> new_data(num_bins_ * new_stride, 0);
        for (size_t bin = 0; bin < num_bins_; ++bin) {
            for (size_t bc = 0; bc < num_barcodes_; ++bc) {
                new_data[bin * new_stride + bc] = data_[bin * stride_ + bc];
            }
        }
        data_ = std::move(new_data);
        stride_ = new_stride;
    }
    num_barcodes_ = n;
}

void CountsMatrix::increment(size_t bin_idx, size_t barcode_idx) {
    if (bin_idx >= num_bins_) {
        throw std::runtime_error("counts matrix: bin index " + std::to_string(bin_idx) +
                                 " out of bounds (num_bins=" + std::to_string(num_bins_) + ")");
    }
    if (barcode_idx >= num_barcodes_) grow_to(barcode_idx + 1);
    data_[bin_idx * stride_ + barcode_idx]++;
}

int32_t CountsMatrix::at(size_t bin_idx, size_t barcode_idx) const {
    if (bin_idx >= num_bins_) {
        throw std::runtime_error("counts matrix: bin index " + std::to_string(bin_idx) +
                                 " out of bounds (num_bins=" + std::to_string(num_bins_) + ")");
    }
    if (barcode_idx >= num_barcodes_) {
        throw std::runtime_error("counts matrix: barcode index " + std::to_string(barcode_idx) +
                                 " out of bounds (num_barcodes=" + std::to_string(num_barcodes_) + ")");
    }
    return data_[bin_idx * stride_ + barcode_idx];
}

void CountsMatrix::write_gz(const std::string& path,
                             const std::vector<std::string>& bin_names,
                             const std::vector<std::string>& barcode_names) const {
    if (bin_names.size() != num_bins_) {
        throw std::runtime_error("counts matrix: bin_names size " + std::to_string(bin_names.size()) +
                                 " does not match num_bins " + std::to_string(num_bins_));
    }
    if (barcode_names.size() != num_barcodes_) {
        throw std::runtime_error("counts matrix: barcode_names size " + std::to_string(barcode_names.size()) +
                                 " does not match num_barcodes " + std::to_string(num_barcodes_));
    }

    gzFile out = gzopen(path.c_str(), "wb");
    if (!out) {
        throw std::runtime_error("cannot open output file for writing: " + path);
    }

    auto checked_write = [&](const std::string& chunk, const std::string& what) {
        int written = gzwrite(out, chunk.data(), static_cast<unsigned>(chunk.size()));
        if (written <= 0 || static_cast<size_t>(written) != chunk.size()) {
            gzclose(out);
            throw std::runtime_error("counts matrix: failed writing " + what + " to " + path);
        }
    };

    std::ostringstream header;
    header << "bin";
    for (const auto& name : barcode_names) header << '\t' << name;
    header << '\n';
    checked_write(header.str(), "header");

    for (size_t bin = 0; bin < num_bins_; ++bin) {
        std::ostringstream row;
        row << bin_names[bin];
        for (size_t bc = 0; bc < barcode_names.size(); ++bc) {
            row << '\t' << at(bin, bc);
        }
        row << '\n';
        checked_write(row.str(), "bin row " + bin_names[bin]);
    }

    int close_ret = gzclose(out);
    if (close_ret != Z_OK) {
        throw std::runtime_error("counts matrix: gzclose failed (code " +
                                  std::to_string(close_ret) + ") for " + path);
    }
}
