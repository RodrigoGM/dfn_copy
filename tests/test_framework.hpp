#pragma once
#include <cstdio>

inline int g_tests_run = 0;
inline int g_tests_failed = 0;

#define ASSERT_TRUE(cond) do { \
    g_tests_run++; \
    if (!(cond)) { \
        std::printf("FAIL %s:%d: ASSERT_TRUE(%s)\n", __FILE__, __LINE__, #cond); \
        g_tests_failed++; \
    } \
} while (0)

#define ASSERT_EQ(a, b) do { \
    g_tests_run++; \
    auto _a = (a); auto _b = (b); \
    if (!(_a == _b)) { \
        std::printf("FAIL %s:%d: ASSERT_EQ(%s, %s)\n", __FILE__, __LINE__, #a, #b); \
        g_tests_failed++; \
    } \
} while (0)

#define TEST_REPORT() do { \
    std::printf("%d assertions, %d failed\n", g_tests_run, g_tests_failed); \
    return g_tests_failed == 0 ? 0 : 1; \
} while (0)

#include <htslib/sam.h>
#include <string>
#include <stdexcept>

// Builds a minimal synthetic alignment record for unit testing: a single
// `seq_len`M CIGAR op, no indels/clips. Caller must bam_destroy1() it.
inline bam1_t* make_test_record(const std::string& qname, uint16_t flag,
                                 int32_t tid, int64_t pos, uint8_t mapq,
                                 int32_t mtid, int64_t mpos, int64_t isize,
                                 int32_t seq_len) {
    bam1_t* rec = bam_init1();
    uint32_t cigar[1] = { static_cast<uint32_t>(bam_cigar_gen(seq_len, BAM_CMATCH)) };
    std::string seq(seq_len, 'A');
    std::string qual(seq_len, static_cast<char>(30));
    int ret = bam_set1(rec, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                        1, cigar, mtid, mpos, isize,
                        seq.size(), seq.c_str(), qual.c_str(), 32);
    if (ret < 0) {
        bam_destroy1(rec);
        throw std::runtime_error("bam_set1 failed for test record " + qname);
    }
    return rec;
}
