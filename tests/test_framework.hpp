#pragma once
#include <cstdio>
#include <cmath>

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

#define ASSERT_NEAR(a, b, tol) do { \
    g_tests_run++; \
    double _a = static_cast<double>(a); double _b = static_cast<double>(b); \
    if (std::fabs(_a - _b) > (tol)) { \
        std::printf("FAIL %s:%d: ASSERT_NEAR(%s, %s) |%.6f - %.6f| > %.6f\n", \
                     __FILE__, __LINE__, #a, #b, _a, _b, static_cast<double>(tol)); \
        g_tests_failed++; \
    } \
} while (0)

#define TEST_REPORT() do { \
    std::printf("%d assertions, %d failed\n", g_tests_run, g_tests_failed); \
    return g_tests_failed == 0 ? 0 : 1; \
} while (0)
