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
