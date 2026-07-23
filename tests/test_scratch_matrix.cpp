#include "test_framework.hpp"
#include "scratch_matrix.hpp"
#include <cstdio>
#include <thread>
#include <vector>

void test_write_then_read_single_cell_round_trips() {
    const char* path = "/tmp/dfn_cbs_test_scratch1.bin";
    auto m = ScratchMatrix::create(path, 3, 4);
    m.write_cell(1, {1.5, 2.5, 3.5, 4.5});
    auto back = m.read_cell(1);
    ASSERT_EQ(back.size(), static_cast<size_t>(4));
    ASSERT_NEAR(back[0], 1.5, 1e-5);
    ASSERT_NEAR(back[3], 4.5, 1e-5);
    m.close();
    std::remove(path);
}

void test_concurrent_writes_to_disjoint_cells_are_safe() {
    const char* path = "/tmp/dfn_cbs_test_scratch2.bin";
    const size_t num_cells = 20;
    const size_t num_bins = 8;
    auto m = ScratchMatrix::create(path, num_cells, num_bins);

    std::vector<std::thread> threads;
    for (size_t c = 0; c < num_cells; ++c) {
        threads.emplace_back([&m, c, num_bins]() {
            std::vector<double> values(num_bins);
            for (size_t b = 0; b < num_bins; ++b) {
                values[b] = static_cast<double>(c * 100 + b);
            }
            m.write_cell(c, values);
        });
    }
    for (auto& t : threads) t.join();

    for (size_t c = 0; c < num_cells; ++c) {
        auto back = m.read_cell(c);
        for (size_t b = 0; b < num_bins; ++b) {
            ASSERT_NEAR(back[b], static_cast<double>(c * 100 + b), 1e-5);
        }
    }
    m.close();
    std::remove(path);
}

void test_read_all_returns_row_major_by_cell() {
    const char* path = "/tmp/dfn_cbs_test_scratch3.bin";
    auto m = ScratchMatrix::create(path, 2, 3);
    m.write_cell(0, {1.0, 2.0, 3.0});
    m.write_cell(1, {4.0, 5.0, 6.0});
    auto flat = m.read_all();
    ASSERT_EQ(flat.size(), static_cast<size_t>(6));
    ASSERT_NEAR(flat[0 * 3 + 0], 1.0f, 1e-5);
    ASSERT_NEAR(flat[0 * 3 + 2], 3.0f, 1e-5);
    ASSERT_NEAR(flat[1 * 3 + 1], 5.0f, 1e-5);
    m.close();
    std::remove(path);
}

int main() {
    test_write_then_read_single_cell_round_trips();
    test_concurrent_writes_to_disjoint_cells_are_safe();
    test_read_all_returns_row_major_by_cell();
    TEST_REPORT();
}
