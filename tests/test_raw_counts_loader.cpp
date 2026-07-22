#include "test_framework.hpp"
#include "raw_counts_loader.hpp"
#include <zlib.h>
#include <cstdio>
#include <fstream>

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

void test_nonexistent_file_throws() {
    bool threw = false;
    try {
        load_raw_counts("/tmp/nonexistent_dfn_cbs_file_xyz.txt.gz", 1, {"chr1:0:1000"});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

void test_genuinely_empty_file_throws() {
    const char* path = "/tmp/dfn_cbs_test_empty.txt.gz";
    write_gz(path, "");

    bool threw = false;
    try {
        load_raw_counts(path, 0, {});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_invalid_header_throws() {
    const char* path = "/tmp/dfn_cbs_test_badheader.txt.gz";
    write_gz(path, "notbin\tAAAA-1\n");

    bool threw = false;
    try {
        load_raw_counts(path, 0, {});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_too_many_bin_rows_throws() {
    const char* path = "/tmp/dfn_cbs_test_toomanyrows.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\n"
        "chr1:0:1000\t10\n"
        "chr1:1000:2000\t20\n");

    bool threw = false;
    try {
        load_raw_counts(path, 1, {"chr1:0:1000"}); // expects 1, file has 2
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_wrong_field_count_throws() {
    const char* path = "/tmp/dfn_cbs_test_wrongcols.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\tCCCC-1\n"
        "chr1:0:1000\t10\n");  // only 1 value, expects 2

    bool threw = false;
    try {
        load_raw_counts(path, 1, {"chr1:0:1000"});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_non_integer_value_throws() {
    const char* path = "/tmp/dfn_cbs_test_nonint.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\n"
        "chr1:0:1000\tnotanumber\n");

    bool threw = false;
    try {
        load_raw_counts(path, 1, {"chr1:0:1000"});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_out_of_range_count_throws() {
    const char* path = "/tmp/dfn_cbs_test_outofrange.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\n"
        "chr1:0:1000\t99999999999\n");  // outside int32 range

    bool threw = false;
    try {
        load_raw_counts(path, 1, {"chr1:0:1000"});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_bin_names_size_mismatch_throws() {
    const char* path = "/tmp/dfn_cbs_test_binnames_mismatch.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\n"
        "chr1:0:1000\t10\n");

    bool threw = false;
    try {
        load_raw_counts(path, 1, {}); // expected_num_bins=1 but empty expected_bin_names
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_too_few_bin_rows_with_matching_bin_names_throws() {
    const char* path = "/tmp/dfn_cbs_test_toofewrows_matching.txt.gz";
    write_gz(path,
        "bin\tAAAA-1\n"
        "chr1:0:1000\t10\n");

    bool threw = false;
    try {
        // Properly-sized expected_bin_names (3 names for 3 bins), but file has only 1 row
        load_raw_counts(path, 3, {"chr1:0:1000", "chr1:1000:2000", "chr1:2000:3000"});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

void test_corrupted_gzip_file_throws() {
    const char* path = "/tmp/dfn_cbs_test_corrupted.txt.gz";

    // Write a valid gzip file first
    write_gz(path, "bin\tAAAA-1\nchr1:0:1000\t10\n");

    // Now truncate it to make it corrupted by writing just the gzip magic bytes
    std::ofstream truncate(path, std::ios::binary | std::ios::out | std::ios::trunc);
    truncate.put(0x1F);
    truncate.put(0x8B);
    truncate.close();

    bool threw = false;
    try {
        load_raw_counts(path, 1, {"chr1:0:1000"});
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        // Verify the error message is about corrupted file, not just "is empty"
        if (msg.find("corrupted") != std::string::npos ||
            msg.find("valid gzip") != std::string::npos) {
            threw = true;
        }
    }
    ASSERT_TRUE(threw);

    std::remove(path);
}

int main() {
    test_loads_matrix_matching_bin_names();
    test_bin_count_mismatch_throws();
    test_bin_name_mismatch_throws();
    test_nonexistent_file_throws();
    test_genuinely_empty_file_throws();
    test_invalid_header_throws();
    test_too_many_bin_rows_throws();
    test_wrong_field_count_throws();
    test_non_integer_value_throws();
    test_out_of_range_count_throws();
    test_bin_names_size_mismatch_throws();
    test_too_few_bin_rows_with_matching_bin_names_throws();
    test_corrupted_gzip_file_throws();
    TEST_REPORT();
}
