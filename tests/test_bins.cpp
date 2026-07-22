#include "test_framework.hpp"
#include "bins.hpp"
#include <fstream>
#include <cstdio>

std::string write_temp_file(const std::string& name, const std::string& content) {
    std::string path = "/tmp/" + name;
    std::ofstream out(path);
    out << content;
    out.close();
    return path;
}

void test_load_basic() {
    std::string path = write_temp_file("bins_basic.tsv",
        "chrom\tstart\tend\tgc\n"
        "chr1\t0\t1000\t0.4\n"
        "chr1\t1000\t2000\t0.5\n"
        "chr2\t0\t1000\t0.45\n");
    BinIndex idx = BinIndex::load(path);
    ASSERT_EQ(idx.size(), static_cast<size_t>(3));
    ASSERT_EQ(idx.name(0), std::string("chr1:0:1000"));
    ASSERT_EQ(idx.name(2), std::string("chr2:0:1000"));
}

void test_load_columns_out_of_order() {
    std::string path = write_temp_file("bins_reorder.tsv",
        "gc\tend\tstart\tchrom\n"
        "0.4\t1000\t0\tchr1\n");
    BinIndex idx = BinIndex::load(path);
    ASSERT_EQ(idx.size(), static_cast<size_t>(1));
    ASSERT_EQ(idx.name(0), std::string("chr1:0:1000"));
}

void test_load_missing_column_throws() {
    std::string path = write_temp_file("bins_missing_col.tsv",
        "chrom\tstart\tend\n"
        "chr1\t0\t1000\n");
    bool threw = false;
    try {
        BinIndex::load(path);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

void test_load_malformed_number_throws() {
    std::string path = write_temp_file("bins_malformed.tsv",
        "chrom\tstart\tend\tgc\n"
        "chr1\tNOTANUMBER\t1000\t0.4\n");
    bool threw = false;
    try {
        BinIndex::load(path);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

void test_find_lookup() {
    std::string path = write_temp_file("bins_lookup.tsv",
        "chrom\tstart\tend\tgc\n"
        "chr1\t0\t1000\t0.4\n"
        "chr1\t1000\t2000\t0.5\n"
        "chr2\t0\t1000\t0.45\n");
    BinIndex idx = BinIndex::load(path);

    auto r1 = idx.find("chr1", 500);
    ASSERT_TRUE(r1.has_value());
    ASSERT_EQ(*r1, static_cast<size_t>(0));

    auto r2 = idx.find("chr1", 1500);
    ASSERT_TRUE(r2.has_value());
    ASSERT_EQ(*r2, static_cast<size_t>(1));

    auto r3 = idx.find("chr1", 999999);
    ASSERT_TRUE(!r3.has_value());

    auto r4 = idx.find("chrX", 10);
    ASSERT_TRUE(!r4.has_value());
}

void test_compute_bin_coordinate() {
    ASSERT_EQ(compute_bin_coordinate(100, 150, PositionMode::Start), static_cast<int64_t>(100));
    ASSERT_EQ(compute_bin_coordinate(100, 150, PositionMode::End), static_cast<int64_t>(149));
    ASSERT_EQ(compute_bin_coordinate(100, 150, PositionMode::Midpoint), static_cast<int64_t>(125));
}

void test_validate_chrom_naming_ok() {
    std::string err = validate_chrom_naming({"chr1", "chr2", "chr3"}, {"chr1", "chr2"});
    ASSERT_EQ(err, std::string(""));
}

void test_validate_chrom_naming_chr_prefix_mismatch() {
    std::string err = validate_chrom_naming({"chr1", "chr2", "chr3"}, {"1", "2"});
    ASSERT_TRUE(!err.empty());
}

void test_validate_chrom_naming_reverse_mismatch() {
    std::string err = validate_chrom_naming({"1", "2", "3"}, {"chr1", "chr2"});
    ASSERT_TRUE(!err.empty());
}

void test_validate_chrom_naming_unknown_chrom() {
    std::string err = validate_chrom_naming({"chr1", "chr2"}, {"chr1", "chrZ"});
    ASSERT_TRUE(!err.empty());
}

int main() {
    test_load_basic();
    test_load_columns_out_of_order();
    test_load_missing_column_throws();
    test_load_malformed_number_throws();
    test_find_lookup();
    test_compute_bin_coordinate();
    test_validate_chrom_naming_ok();
    test_validate_chrom_naming_chr_prefix_mismatch();
    test_validate_chrom_naming_reverse_mismatch();
    test_validate_chrom_naming_unknown_chrom();
    TEST_REPORT();
}
