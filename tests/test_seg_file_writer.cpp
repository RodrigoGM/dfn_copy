#include "test_framework.hpp"
#include "seg_file_writer.hpp"
#include <fstream>
#include <sstream>
#include <thread>
#include <cstdio>
#include <set>

void test_writes_header_and_rows() {
    const char* path = "/tmp/dfn_cbs_test.seg";
    {
        SegFileWriter w(path);
        w.append_rows({{"AAAA-1", "chr1", 0, 1000, 1, 1.0}});
    }
    std::ifstream f(path);
    std::string header;
    std::getline(f, header);
    ASSERT_EQ(header, std::string("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean"));
    std::string row;
    std::getline(f, row);
    ASSERT_EQ(row, std::string("AAAA-1\tchr1\t0\t1000\t1\t1"));
    std::remove(path);
}

void test_concurrent_appends_are_not_corrupted() {
    const char* path = "/tmp/dfn_cbs_test_concurrent.seg";
    {
        SegFileWriter w(path);
        std::vector<std::thread> threads;
        for (int t = 0; t < 10; ++t) {
            threads.emplace_back([&w, t]() {
                std::string barcode = "CELL" + std::to_string(t);
                w.append_rows({{barcode, "chr1", 0, 1000, 1, static_cast<double>(t)}});
            });
        }
        for (auto& th : threads) th.join();
    }

    std::ifstream f(path);
    std::string line;
    std::getline(f, line); // header
    std::set<std::string> seen_barcodes;
    int row_count = 0;
    while (std::getline(f, line)) {
        ++row_count;
        std::istringstream iss(line);
        std::string barcode;
        std::getline(iss, barcode, '\t');
        seen_barcodes.insert(barcode);
    }
    ASSERT_EQ(row_count, 10);
    ASSERT_EQ(seen_barcodes.size(), static_cast<size_t>(10));
    std::remove(path);
}

void test_unwritable_path_throws() {
    bool threw = false;
    try {
        SegFileWriter w("/tmp/dfn_cbs_nonexistent_dir/out.seg");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    ASSERT_TRUE(threw);
}

int main() {
    test_writes_header_and_rows();
    test_concurrent_appends_are_not_corrupted();
    test_unwritable_path_throws();
    TEST_REPORT();
}
