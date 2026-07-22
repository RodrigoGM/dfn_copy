#include "test_framework.hpp"
#include "fragment_pairing.hpp"

FragmentCandidate make_candidate(std::string qname, int32_t tid, int64_t start,
                                  int64_t end, int64_t mate_pos, size_t barcode_idx,
                                  bool is_read1) {
    FragmentCandidate c;
    c.qname = qname;
    c.tid = tid;
    c.align_start = start;
    c.align_end = end;
    c.mate_pos = mate_pos;
    c.barcode_idx = barcode_idx;
    c.is_read1 = is_read1;
    return c;
}

void test_auto_prefers_longer_mate() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);   // 27bp
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);  // 51bp

    auto result1 = buf.process(r1);
    ASSERT_TRUE(!result1.has_value());
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(1));

    auto result2 = buf.process(r2);
    ASSERT_TRUE(result2.has_value());
    ASSERT_EQ(result2->align_start, static_cast<int64_t>(520));
    ASSERT_TRUE(!result2->is_read1);
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(0));
}

void test_auto_tie_breaks_to_read1() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    auto r2 = make_candidate("q1", 0, 500, 550, 600, 0, false);
    auto r1 = make_candidate("q1", 0, 600, 650, 500, 0, true);

    buf.process(r2);
    auto result = buf.process(r1);
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result->is_read1);
}

void test_r1_mode_picks_fixed_mate() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R1);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);

    buf.process(r1);
    auto result = buf.process(r2);
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result->is_read1);
    ASSERT_EQ(result->align_start, static_cast<int64_t>(500));
}

void test_r2_mode_picks_fixed_mate_when_both_present() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R2);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);

    buf.process(r1);
    auto result = buf.process(r2);
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(!result->is_read1);
}

void test_auto_singleton_flush_is_counted() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    buf.process(r1);

    auto flushed = buf.flush_up_to(0, 300); // hasn't passed mate_pos=520 yet
    ASSERT_EQ(flushed.size(), static_cast<size_t>(0));
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(1));

    auto flushed2 = buf.flush_up_to(0, 600); // now past mate_pos=520
    ASSERT_EQ(flushed2.size(), static_cast<size_t>(1));
    ASSERT_EQ(flushed2[0].align_start, static_cast<int64_t>(500));
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(0));
}

void test_r1_mode_singleton_read2_only_is_dropped() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R1);
    auto r2 = make_candidate("q1", 0, 520, 571, 500, 0, false);
    buf.process(r2);

    auto flushed = buf.flush_up_to(0, 600);
    ASSERT_EQ(flushed.size(), static_cast<size_t>(0));
}

void test_r1_mode_singleton_read1_only_is_kept() {
    FragmentPairingBuffer buf(SingleEndCountingMode::R1);
    auto r1 = make_candidate("q1", 0, 500, 527, 520, 0, true);
    buf.process(r1);

    auto flushed = buf.flush_up_to(0, 600);
    ASSERT_EQ(flushed.size(), static_cast<size_t>(1));
    ASSERT_TRUE(flushed[0].is_read1);
}

void test_flush_all_at_end_of_stream() {
    FragmentPairingBuffer buf(SingleEndCountingMode::Auto);
    buf.process(make_candidate("q1", 0, 500, 527, 520, 0, true));
    buf.process(make_candidate("q2", 1, 100, 150, 200, 1, true));

    auto flushed = buf.flush_all();
    ASSERT_EQ(flushed.size(), static_cast<size_t>(2));
    ASSERT_EQ(buf.pending_count(), static_cast<size_t>(0));
}

int main() {
    test_auto_prefers_longer_mate();
    test_auto_tie_breaks_to_read1();
    test_r1_mode_picks_fixed_mate();
    test_r2_mode_picks_fixed_mate_when_both_present();
    test_auto_singleton_flush_is_counted();
    test_r1_mode_singleton_read2_only_is_dropped();
    test_r1_mode_singleton_read1_only_is_kept();
    test_flush_all_at_end_of_stream();
    TEST_REPORT();
}
