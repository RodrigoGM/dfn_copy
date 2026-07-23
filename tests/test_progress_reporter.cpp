#include "test_framework.hpp"
#include "progress_reporter.hpp"
#include <thread>
#include <vector>

void test_format_duration() {
    ASSERT_EQ(format_duration(0.0), std::string("0:00:00"));
    ASSERT_EQ(format_duration(59.0), std::string("0:00:59"));
    ASSERT_EQ(format_duration(3661.0), std::string("1:01:01"));
    ASSERT_EQ(format_duration(-5.0), std::string("0:00:00")); // clamped, never negative
}

void test_format_with_commas() {
    ASSERT_EQ(format_with_commas(0), std::string("0"));
    ASSERT_EQ(format_with_commas(999), std::string("999"));
    ASSERT_EQ(format_with_commas(1000), std::string("1,000"));
    ASSERT_EQ(format_with_commas(1234567), std::string("1,234,567"));
}

void test_format_phase2_status_exact_content() {
    // completed=100, total=400, elapsed=10s -> rate=10 cells/s,
    // remaining=300 cells -> ETA = 300/10 = 30s.
    std::string status = format_phase2_status(100, 400, 10.0);
    ASSERT_EQ(status, std::string(
        "Phase 2: 100 / 400 cells (25.0%) | elapsed 0:00:10 | ETA 0:00:30 | 10.0 cells/s"));
}

void test_format_phase2_status_handles_zero_elapsed() {
    // Guards against divide-by-zero on the very first tick.
    std::string status = format_phase2_status(0, 100, 0.0);
    ASSERT_TRUE(status.find("0.0 cells/s") != std::string::npos);
}

void test_reporter_increment_from_multiple_threads() {
    Phase2Reporter reporter(1000, /*quiet=*/true); // quiet: no background thread/printing
    std::vector<std::thread> threads;
    for (int t = 0; t < 10; ++t) {
        threads.emplace_back([&reporter]() {
            for (int i = 0; i < 100; ++i) reporter.increment();
        });
    }
    for (auto& th : threads) th.join();
    ASSERT_EQ(reporter.completed(), static_cast<size_t>(1000));
    reporter.stop(); // idempotent, must not hang
}

int main() {
    test_format_duration();
    test_format_with_commas();
    test_format_phase2_status_exact_content();
    test_format_phase2_status_handles_zero_elapsed();
    test_reporter_increment_from_multiple_threads();
    TEST_REPORT();
}
