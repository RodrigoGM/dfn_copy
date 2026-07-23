#include "test_framework.hpp"
#include "transpose_writer.hpp"
#include <zlib.h>
#include <cstdio>
#include <sstream>
#include <csignal>
#include <sys/resource.h>
#include <stdexcept>

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

// Real failure-injection test: cap RLIMIT_FSIZE so the underlying write()
// syscall gzwrite() ultimately performs fails partway through with EFBIG,
// exercising the checked_write()/gzclose() error paths against an actual
// short/failed write rather than a mock. SIGXFSZ is ignored so the process
// gets the EFBIG failure back through the return value instead of being
// killed. RLIMIT_FSIZE is a process-wide resource limit, restored in a RAII
// guard so it never leaks into other tests in this binary.
struct FsizeLimitGuard {
    struct rlimit old_lim{};
    struct sigaction old_sigxfsz{};
    FsizeLimitGuard(rlim_t cur_limit) {
        getrlimit(RLIMIT_FSIZE, &old_lim);
        struct sigaction ignore_action{};
        ignore_action.sa_handler = SIG_IGN;
        sigaction(SIGXFSZ, &ignore_action, &old_sigxfsz);
        struct rlimit lim = old_lim;
        lim.rlim_cur = cur_limit;
        if (setrlimit(RLIMIT_FSIZE, &lim) != 0) {
            throw std::runtime_error("test setup: setrlimit(RLIMIT_FSIZE) failed");
        }
    }
    ~FsizeLimitGuard() {
        setrlimit(RLIMIT_FSIZE, &old_lim);
        sigaction(SIGXFSZ, &old_sigxfsz, nullptr);
    }
};

void test_gzwrite_failure_mid_stream_throws() {
    // Enough distinct, high-entropy float values that gzip can't compress
    // its way under the byte cap -- forces the cap to bind partway through
    // the row loop rather than never being reached.
    const size_t num_cells = 300;
    const size_t num_bins = 50;
    std::vector<float> scratch_data(num_cells * num_bins);
    for (size_t i = 0; i < scratch_data.size(); ++i) {
        // Simple LCG-derived spread, not just a repeating pattern, so the
        // deflate stream stays incompressible enough to hit the cap.
        scratch_data[i] = static_cast<float>((i * 2654435761u) % 1000003u) / 7.0f;
    }
    std::vector<std::string> bin_names;
    for (size_t b = 0; b < num_bins; ++b) bin_names.push_back("chr1:" + std::to_string(b * 1000) + ":" + std::to_string(b * 1000 + 1000));
    std::vector<std::string> barcode_names;
    for (size_t c = 0; c < num_cells; ++c) barcode_names.push_back("BARCODE" + std::to_string(c) + "-1");

    const char* path = "/tmp/dfn_cbs_test_transpose_rlimit.txt.gz";
    std::remove(path);

    bool threw = false;
    std::string what;
    {
        FsizeLimitGuard guard(200); // bytes -- far smaller than the full output
        try {
            write_transposed_gz(path, scratch_data, num_cells, num_bins, bin_names, barcode_names);
        } catch (const std::runtime_error& e) {
            threw = true;
            what = e.what();
        }
    }
    ASSERT_TRUE(threw);
    // Confirm it's actually our checked-write/gzclose error, not some
    // unrelated failure, and that it mentions the path.
    ASSERT_TRUE(what.find("transpose writer") != std::string::npos);
    ASSERT_TRUE(what.find(path) != std::string::npos);
    std::remove(path);
}

int main() {
    test_writes_bins_by_cells_layout();
    test_size_mismatch_throws();
    test_gzwrite_failure_mid_stream_throws();
    TEST_REPORT();
}
