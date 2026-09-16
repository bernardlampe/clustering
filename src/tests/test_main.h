#ifndef __TEST_MAIN_H__
#define __TEST_MAIN_H__

/* Minimal assertion-based test framework shared by the test binaries.
 * CHECK(expr)              - boolean assertion
 * CHECK_NEAR(a, b, eps)    - floating point comparison
 * CHECK_THROWS(expr)       - expression must throw (Exception or any type)
 * test_report()            - print summary, return non-zero on failure
 */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "Exception.h"

static int g_pass = 0;
static int g_fail = 0;
static std::vector<std::string> g_failures;

static inline void check_(bool ok, const char *expr, const char *file, int line) {
  if (ok) {
    g_pass++;
  } else {
    g_fail++;
    char buf[1024];
    snprintf(buf, sizeof buf, "%s:%d: %s", file, line, expr);
    g_failures.push_back(buf);
  }
}

#define CHECK(expr) check_((expr), #expr, __FILE__, __LINE__)

// yields bool so it can compose inside CHECK(...) expressions
#define CHECK_NEAR(a, b, eps)                                                  \
  test_near_impl_((a), (b), (eps), #a " == " #b, __FILE__, __LINE__)

static inline bool test_near_impl_(double a, double b, double eps, const char *expr,
                                   const char *file, int line) {
  if (std::fabs(a - b) <= eps) {
    g_pass++;
    return true;
  }
  g_fail++;
  char buf[1024];
  snprintf(buf, sizeof buf, "%s:%d: %s (got %.8g want %.8g)", file, line, expr, a, b);
  g_failures.push_back(buf);
  return false;
}
#define CHECK_THROWS(expr)                                                     \
  do {                                                                         \
    bool threw_ = false;                                                       \
    try {                                                                      \
      expr;                                                                    \
    } catch (Exception &) {                                                    \
      threw_ = true;                                                           \
    } catch (...) {                                                            \
      threw_ = true;                                                           \
    }                                                                          \
    check_(threw_, #expr " must throw", __FILE__, __LINE__);                   \
  } while (0)

static inline int test_report() {
  for (size_t i = 0; i < g_failures.size(); i++)
    fprintf(stderr, "FAIL %s\n", g_failures[i].c_str());
  fprintf(stderr, "%d passed, %d failed\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}

#endif // __TEST_MAIN_H__
