#include "test_framework.hpp"
#include "transpose_writer.hpp"
#include <zlib.h>
#include <cstdio>
#include <sstream>

std::string read_gz(const std::string& path) {
    gzFile in = gzopen(path.c_str(), "rb");
    std::string content;
    char buf[256];
    int n;
    while ((n = gzread(in, buf, sizeof(buf))) > 0) content.append(buf, n);
    gzclose(in);
    return content;
}

void test_writes_bins_by_cells_layout() {
    // 2 cells x 3 bins, cell-major (row=cell): cell0={1,2,3}, cell1={4,5,6}.
    std::vector<float> scratch_data = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    const char* path = "/tmp/dfn_cbs_test_transpose.txt.gz";
    write_transposed_gz(path, scratch_data, 2, 3,
                         {"chr1:0:1000", "chr1:1000:2000", "chr2:0:1000"},
                         {"AAAA-1", "CCCC-1"});

    std::string content = read_gz(path);
    std::istringstream iss(content);
    std::string line;
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("bin\tAAAA-1\tCCCC-1"));
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("chr1:0:1000\t1\t4"));
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("chr1:1000:2000\t2\t5"));
    std::getline(iss, line);
    ASSERT_EQ(line, std::string("chr2:0:1000\t3\t6"));

    std::remove(path);
}

void test_size_mismatch_throws() {
    std::vector<float> scratch_data = {1.0f, 2.0f}; // wrong size for 2x3
    bool threw = false;
    try {
        write_transposed_gz("/tmp/dfn_cbs_test_transpose_bad.txt.gz", scratch_data, 2, 3,
                             {"a", "b", "c"}, {"X", "Y"});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

int main() {
    test_writes_bins_by_cells_layout();
    test_size_mismatch_throws();
    TEST_REPORT();
}
