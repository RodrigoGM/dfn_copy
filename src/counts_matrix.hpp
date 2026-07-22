#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

class CountsMatrix {
public:
    explicit CountsMatrix(size_t num_bins, size_t initial_barcodes = 0);

    void grow_to(size_t num_barcodes);
    void increment(size_t bin_idx, size_t barcode_idx);

    void write_gz(const std::string& path,
                  const std::vector<std::string>& bin_names,
                  const std::vector<std::string>& barcode_names) const;

    int32_t at(size_t bin_idx, size_t barcode_idx) const;
    size_t num_bins() const { return num_bins_; }
    size_t num_barcodes() const { return num_barcodes_; }

private:
    size_t num_bins_;
    size_t num_barcodes_;
    size_t stride_;
    std::vector<int32_t> data_;
};
