#include "fragment_pairing.hpp"

FragmentPairingBuffer::FragmentPairingBuffer(SingleEndCountingMode mode) : mode_(mode) {}

std::optional<FragmentCandidate> FragmentPairingBuffer::process(
        const FragmentCandidate& incoming) {
    auto it = pending_.find(incoming.qname);
    if (it == pending_.end()) {
        pending_.emplace(incoming.qname, incoming);
        return std::nullopt;
    }
    FragmentCandidate other = it->second;
    pending_.erase(it);
    return resolve_pair(incoming, other);
}

std::optional<FragmentCandidate> FragmentPairingBuffer::resolve_pair(
        const FragmentCandidate& a, const FragmentCandidate& b) const {
    switch (mode_) {
        case SingleEndCountingMode::Auto: {
            int64_t len_a = a.align_end - a.align_start;
            int64_t len_b = b.align_end - b.align_start;
            if (len_a != len_b) return (len_a > len_b) ? a : b;
            return a.is_read1 ? a : b;
        }
        case SingleEndCountingMode::R1:
            if (a.is_read1) return a;
            if (b.is_read1) return b;
            return std::nullopt;
        case SingleEndCountingMode::R2:
            if (!a.is_read1) return a;
            if (!b.is_read1) return b;
            return std::nullopt;
    }
    return std::nullopt;
}

std::optional<FragmentCandidate> FragmentPairingBuffer::resolve_singleton(
        const FragmentCandidate& c) const {
    switch (mode_) {
        case SingleEndCountingMode::Auto:
            return c;
        case SingleEndCountingMode::R1:
            return c.is_read1 ? std::optional<FragmentCandidate>(c) : std::nullopt;
        case SingleEndCountingMode::R2:
            return !c.is_read1 ? std::optional<FragmentCandidate>(c) : std::nullopt;
    }
    return std::nullopt;
}

std::vector<FragmentCandidate> FragmentPairingBuffer::flush_up_to(
        int32_t tid, int64_t before_pos) {
    std::vector<FragmentCandidate> results;
    for (auto it = pending_.begin(); it != pending_.end(); ) {
        const FragmentCandidate& c = it->second;
        if (c.tid == tid && c.mate_pos < before_pos) {
            auto resolved = resolve_singleton(c);
            if (resolved) results.push_back(*resolved);
            it = pending_.erase(it);
        } else {
            ++it;
        }
    }
    return results;
}

std::vector<FragmentCandidate> FragmentPairingBuffer::flush_all() {
    std::vector<FragmentCandidate> results;
    for (auto& entry : pending_) {
        auto resolved = resolve_singleton(entry.second);
        if (resolved) results.push_back(*resolved);
    }
    pending_.clear();
    return results;
}
