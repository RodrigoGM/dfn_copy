#include "test_framework.hpp"
#include "counts_matrix.hpp"
#include <zlib.h>
#include <sstream>

void test_increment_and_read_back() {
    CountsMatrix m(2, 2);
    m.increment(0, 0);
    m.increment(0, 0);
    m.increment(1, 1);
    ASSERT_EQ(m.at(0, 0), 2);
    ASSERT_EQ(m.at(0, 1), 0);
    ASSERT_EQ(m.at(1, 1), 1);
}

void test_grow_preserves_existing_data() {
    CountsMatrix m(2, 1);
    m.increment(0, 0);
    m.increment(1, 0);
    m.grow_to(3);
    ASSERT_EQ(m.num_barcodes(), static_cast<size_t>(3));
    ASSERT_EQ(m.at(0, 0), 1);
    ASSERT_EQ(m.at(1, 0), 1);
    ASSERT_EQ(m.at(0, 2), 0);
    m.increment(0, 2);
    ASSERT_EQ(m.at(0, 2), 1);
}

void test_increment_grows_implicitly() {
    CountsMatrix m(1, 0);
    m.increment(0, 5);
    ASSERT_EQ(m.num_barcodes(), static_cast<size_t>(6));
    ASSERT_EQ(m.at(0, 5), 1);
}

std::string read_gz_file(const std::string& path) {
    gzFile f = gzopen(path.c_str(), "rb");
    std::string content;
    char buf[4096];
    int n;
    while ((n = gzread(f, buf, sizeof(buf))) > 0) {
        content.append(buf, n);
    }
    gzclose(f);
    return content;
}

void test_write_gz_format() {
    CountsMatrix m(2, 2);
    m.increment(0, 0);
    m.increment(0, 0);
    m.increment(1, 1);
    std::string path = "/tmp/test_counts_matrix_out.txt.gz";
    m.write_gz(path, {"chr1:0:1000", "chr1:1000:2000"}, {"AAAA-1", "CCCC-1"});

    std::string content = read_gz_file(path);
    std::string expected =
        "bin\tAAAA-1\tCCCC-1\n"
        "chr1:0:1000\t2\t0\n"
        "chr1:1000:2000\t0\t1\n";
    ASSERT_EQ(content, expected);
}

int main() {
    test_increment_and_read_back();
    test_grow_preserves_existing_data();
    test_increment_grows_implicitly();
    test_write_gz_format();
    TEST_REPORT();
}
