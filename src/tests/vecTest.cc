#include <cmath>
#include <iostream>

#include "test_main.h"
#include "Vec.h"

static void testArithmetic() {
  float a[] = {1, 2, 3};
  float b[] = {4, 5, 6};
  Vec<float> v(a, 3), w(b, 3);

  Vec<float> s = v + w;
  CHECK(s[0] == 5 && s[1] == 7 && s[2] == 9);
  Vec<float> d = w - v;
  CHECK(d[0] == 3 && d[1] == 3 && d[2] == 3);
  Vec<float> p = v * w;
  CHECK(p[0] == 4 && p[1] == 10 && p[2] == 18);
  Vec<float> q = w / v;
  CHECK(q[0] == 4 && std::fabs(q[1] - 2.5f) < 1e-6f && q[2] == 2);

  Vec<float> r = v;
  r += w; CHECK(r[0] == 5 && r[1] == 7 && r[2] == 9);
  r = v; r *= w; CHECK(r[0] == 4 && r[1] == 10 && r[2] == 18);
  r = w; r /= v; CHECK(r[0] == 4 && CHECK_NEAR(r[1], 2.5, 1e-6) && r[2] == 2);
}

static void testScalarOps() {
  float a[] = {1, 2, 3};
  Vec<float> v(a, 3);
  Vec<float> s = v * 2.0f;
  CHECK(s[0] == 2 && s[1] == 4 && s[2] == 6);
  Vec<float> t = v + 1.0f;
  CHECK(t[0] == 2 && t[1] == 3 && t[2] == 4);
  Vec<float> u = v - 1.0f;
  CHECK(u[0] == 0 && u[1] == 1 && u[2] == 2);
  Vec<float> x = v / 2.0f;
  CHECK_NEAR(x[0], 0.5, 1e-7);
  v += 1.0f; v -= 1.0f; v *= 2.0f; v /= 2.0f;
  CHECK(v[0] == 1 && v[1] == 2 && v[2] == 3);
}

static void testReductions() {
  float a[] = {1, 2, 3, 4};
  Vec<float> v(a, 4);
  CHECK_NEAR(v.sum(), 10.0, 1e-6);
  CHECK_NEAR(v.prod(), 24.0, 1e-6);  // regression: prod() was seeded with 0
  CHECK_NEAR(v.norm(), std::sqrt(30.0), 1e-6);
  float b[] = {-1, 2, -3};
  Vec<float> w(b, 3);
  Vec<float> aw = w.abs();
  CHECK(aw[0] == 1 && aw[1] == 2 && aw[2] == 3);
  Vec<float> z; // empty
  CHECK(z.prod() == 0 && z.sum() == 0);
}

static void testDimensionMismatchThrows() {
  // regression: these threw raw const char*, uncatchable as Exception
  float a[] = {1, 2, 3};
  Vec<float> v(a, 3), w(a, 2);
  CHECK_THROWS(v + w);
  CHECK_THROWS(v - w);
  CHECK_THROWS(v * w);
  CHECK_THROWS(v / w);
  CHECK_THROWS(v += w);
  CHECK_THROWS(v -= w);
  CHECK_THROWS(v *= w);
  CHECK_THROWS(v /= w);
}

int main() {
  testArithmetic();
  testScalarOps();
  testReductions();
  testDimensionMismatchThrows();
  return test_report();
}
