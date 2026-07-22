#pragma once
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include <optional>
#include "cli_args.hpp"

struct FragmentCandidate {
    std::string qname;
    int32_t tid = -1;
    int64_t align_start = 0;
    int64_t align_end = 0;
    int64_t mate_pos = -1;
    size_t barcode_idx = 0;
    bool is_read1 = true;
};

class FragmentPairingBuffer {
public:
    explicit FragmentPairingBuffer(SingleEndCountingMode mode);

    std::optional<FragmentCandidate> process(const FragmentCandidate& incoming);
    std::vector<FragmentCandidate> flush_up_to(int32_t tid, int64_t before_pos);
    std::vector<FragmentCandidate> flush_all();

    size_t pending_count() const { return pending_.size(); }

private:
    SingleEndCountingMode mode_;
    std::unordered_map<std::string, FragmentCandidate> pending_;

    std::optional<FragmentCandidate> resolve_pair(const FragmentCandidate& a,
                                                   const FragmentCandidate& b) const;
    std::optional<FragmentCandidate> resolve_singleton(const FragmentCandidate& c) const;
};
