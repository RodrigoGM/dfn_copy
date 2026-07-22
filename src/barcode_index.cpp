#include "barcode_index.hpp"
#include <fstream>
#include <stdexcept>

BarcodeIndex BarcodeIndex::load_allowlist(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot read barcode allowlist file: " + path);
    }

    BarcodeIndex idx;
    idx.has_allowlist_ = true;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (idx.index_.find(line) != idx.index_.end()) continue;
        idx.index_[line] = idx.names_.size();
        idx.names_.push_back(line);
    }

    if (idx.names_.empty()) {
        throw std::runtime_error("barcode allowlist file is empty: " + path);
    }

    return idx;
}

std::optional<size_t> BarcodeIndex::get_or_create(const std::string& barcode) {
    auto it = index_.find(barcode);
    if (it != index_.end()) return it->second;

    if (has_allowlist_) return std::nullopt;

    size_t new_index = names_.size();
    index_[barcode] = new_index;
    names_.push_back(barcode);
    return new_index;
}
