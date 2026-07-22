#include "test_framework.hpp"
#include "raw_counts_loader.hpp"
#include <zlib.h>
#include <cstdio>

namespace {
void write_gz(const std::string& path, const std::string& content) {
    gzFile out = gzopen(path.c_str(), "wb");
    gzwrite(out, content.data(), static_cast<unsigned>(content.size()));
    gzclose(out);
}
}

void test_loads_matrix_matching_bin_names() {
    const char* path = "/tmp/dfn_cbs_test_counts.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\tCCCC-1\n"
        "chr1:0:1000\t10\t20\n"
        "chr1:1000:2000\t5\t0\n");

    auto m = load_raw_counts(path, 2, {"chr1:0:1000", "chr1:1000:2000"});
    ASSERT_EQ(m.num_bins, static_cast<size_t>(2));
    ASSERT_EQ(m.num_cells, static_cast<size_t>(2));
    ASSERT_EQ(m.barcodes[0], std::string("AAAA-1"));
    ASSERT_EQ(m.barcodes[1], std::string("CCCC-1"));
    ASSERT_EQ(m.data[0 * m.num_cells + 0], static_cast<int32_t>(10));
    ASSERT_EQ(m.data[0 * m.num_cells + 1], static_cast<int32_t>(20));
    ASSERT_EQ(m.data[1 * m.num_cells + 0], static_cast<int32_t>(5));
    ASSERT_EQ(m.data[1 * m.num_cells + 1], static_cast<int32_t>(0));

    std::remove(path);
}

void test_bin_count_mismatch_throws() {
    const char* path = "/tmp/dfn_cbs_test_counts_mismatch.txt.gz";
    write_gz(path, "bin\tAAAA-1\nchr1:0:1000\t10\n");

    bool threw = false;
    try {
        load_raw_counts(path, 5, {"chr1:0:1000"}); // expects 5 bins, file has 1
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_bin_name_mismatch_throws() {
    const char* path = "/tmp/dfn_cbs_test_counts_namemismatch.txt.gz";
    write_gz(path, "bin\tAAAA-1\nchr1:0:1000\t10\n");

    bool threw = false;
    try {
        load_raw_counts(path, 1, {"chr2:0:1000"}); // right count, wrong name/order
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

int main() {
    test_loads_matrix_matching_bin_names();
    test_bin_count_mismatch_throws();
    test_bin_name_mismatch_throws();
    TEST_REPORT();
}
