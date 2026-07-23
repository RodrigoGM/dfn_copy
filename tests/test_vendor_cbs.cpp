#include "test_framework.hpp"
#include "cbs.hpp"
#include <algorithm>

void test_detects_an_injected_breakpoint() {
    // 60 points: first 30 at ~1.0, last 30 at ~2.0 (a clean step). With
    // centering OFF (dfn_cbs always calls it this way), the segment means
    // must reflect the true absolute levels, not a chromosome-recentered
    // difference.
    Series S;
    S.chrom = "chr1";
    for (int i = 0; i < 60; ++i) {
        S.pos.push_back(i * 1000);
        S.x.push_back(i < 30 ? 1.0 : 2.0);
    }

    Args a;
    a.alpha = 0.05;
    a.perms = 200;
    a.min_seg_len = 5;
    a.max_depth = 10;
    a.method = "1cp";
    a.center = false; // dfn_cbs always forces this off

    std::vector<Segment> segs;
    std::mt19937_64 rng(12345);
    segment_chromosome(S, a, segs, rng);
    std::sort(segs.begin(), segs.end(), [](const Segment& x, const Segment& y) { return x.s < y.s; });

    ASSERT_TRUE(segs.size() >= 2);
    // First segment should end near index 29 (the injected breakpoint).
    ASSERT_TRUE(segs.front().e >= 25 && segs.front().e <= 34);
    // Segment means should reflect the true absolute levels (1.0 / 2.0),
    // not values shifted toward zero by centering.
    ASSERT_NEAR(segs.front().mean, 1.0, 0.2);
    ASSERT_NEAR(segs.back().mean, 2.0, 0.2);
}

int main() {
    test_detects_an_injected_breakpoint();
    TEST_REPORT();
}
