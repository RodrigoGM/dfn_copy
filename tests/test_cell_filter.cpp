#include "test_framework.hpp"
#include "cell_filter.hpp"

RawCountsMatrix make_matrix(size_t num_bins, size_t num_cells,
                             std::vector<int32_t> data,
                             std::vector<std::string> barcodes) {
    RawCountsMatrix m;
    m.num_bins = num_bins;
    m.num_cells = num_cells;
    m.data = std::move(data);
    m.barcodes = std::move(barcodes);
    return m;
}

void test_drops_cells_below_threshold() {
    // 2 bins x 3 cells. Column sums: cell0=10+5=15, cell1=100+50=150, cell2=0+0=0.
    auto m = make_matrix(2, 3,
        {10, 100, 0,
          5,  50, 0},
        {"AAAA-1", "CCCC-1", "GGGG-1"});

    auto result = filter_cells(m, 20);
    ASSERT_EQ(result.surviving_indices.size(), static_cast<size_t>(1));
    ASSERT_EQ(result.surviving_indices[0], static_cast<size_t>(1)); // only cell1 (150) survives
    ASSERT_EQ(result.max_column_sum, static_cast<int64_t>(150));
}

void test_all_cells_survive_when_threshold_is_zero() {
    auto m = make_matrix(1, 2, {0, 1}, {"AAAA-1", "CCCC-1"});
    auto result = filter_cells(m, 0);
    ASSERT_EQ(result.surviving_indices.size(), static_cast<size_t>(2));
}

void test_no_cells_survive_reports_max_sum() {
    auto m = make_matrix(1, 2, {5, 9}, {"AAAA-1", "CCCC-1"});
    auto result = filter_cells(m, 100000);
    ASSERT_EQ(result.surviving_indices.size(), static_cast<size_t>(0));
    ASSERT_EQ(result.max_column_sum, static_cast<int64_t>(9));
}

int main() {
    test_drops_cells_below_threshold();
    test_all_cells_survive_when_threshold_is_zero();
    test_no_cells_survive_reports_max_sum();
    TEST_REPORT();
}
