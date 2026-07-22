#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <optional>

class BarcodeIndex {
public:
    BarcodeIndex() = default;

    static BarcodeIndex load_allowlist(const std::string& path);

    bool has_allowlist() const { return has_allowlist_; }
    std::optional<size_t> get_or_create(const std::string& barcode);
    size_t size() const { return names_.size(); }
    const std::vector<std::string>& names() const { return names_; }

private:
    bool has_allowlist_ = false;
    std::vector<std::string> names_;
    std::unordered_map<std::string, size_t> index_;
};
