#include "test_framework.hpp"
#include "bin_gc.hpp"
#include <fstream>
#include <cstdio>

void test_loads_bins_in_file_order() {
    const char* path = "/tmp/dfn_cbs_test_bins.tsv";
    std::ofstream f(path);
    f << "chrom\tstart\tend\tgc\n";
    f << "chr1\t0\t1000\t0.42\n";
    f << "chr1\t1000\t2000\t0.55\n";
    f << "chr2\t0\t1000\t0.60\n";
    f.close();

    auto bins = load_bin_gc(path);
    ASSERT_EQ(bins.size(), static_cast<size_t>(3));
    ASSERT_EQ(bins[0].chrom, std::string("chr1"));
    ASSERT_EQ(bins[0].start, static_cast<int64_t>(0));
    ASSERT_EQ(bins[0].end, static_cast<int64_t>(1000));
    ASSERT_NEAR(bins[0].gc, 0.42, 1e-9);
    ASSERT_EQ(bins[2].chrom, std::string("chr2"));
    ASSERT_EQ(bin_name(bins[0]), std::string("chr1:0:1000"));

    std::remove(path);
}

void test_columns_can_be_reordered() {
    const char* path = "/tmp/dfn_cbs_test_bins_reorder.tsv";
    std::ofstream f(path);
    f << "gc\tchrom\tend\tstart\n";
    f << "0.5\tchrX\t2000\t1000\n";
    f.close();

    auto bins = load_bin_gc(path);
    ASSERT_EQ(bins.size(), static_cast<size_t>(1));
    ASSERT_EQ(bins[0].chrom, std::string("chrX"));
    ASSERT_EQ(bins[0].start, static_cast<int64_t>(1000));
    ASSERT_EQ(bins[0].end, static_cast<int64_t>(2000));

    std::remove(path);
}

void test_missing_column_throws() {
    const char* path = "/tmp/dfn_cbs_test_bins_missing.tsv";
    std::ofstream f(path);
    f << "chrom\tstart\tend\n"; // no gc column
    f << "chr1\t0\t1000\n";
    f.close();

    bool threw = false;
    try {
        load_bin_gc(path);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

int main() {
    test_loads_bins_in_file_order();
    test_columns_can_be_reordered();
    test_missing_column_throws();
    TEST_REPORT();
}
