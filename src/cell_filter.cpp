#include "cell_filter.hpp"

FilteredCells filter_cells(const RawCountsMatrix& m, int64_t min_reads) {
    FilteredCells result;
    result.max_column_sum = 0;

    std::vector<int64_t> column_sums(m.num_cells, 0);
    for (size_t bin = 0; bin < m.num_bins; ++bin) {
        for (size_t cell = 0; cell < m.num_cells; ++cell) {
            column_sums[cell] += m.data[bin * m.num_cells + cell];
        }
    }

    for (size_t cell = 0; cell < m.num_cells; ++cell) {
        if (column_sums[cell] > result.max_column_sum) {
            result.max_column_sum = column_sums[cell];
        }
        if (column_sums[cell] >= min_reads) {
            result.surviving_indices.push_back(cell);
        }
    }

    return result;
}
