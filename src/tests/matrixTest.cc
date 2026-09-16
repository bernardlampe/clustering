#include <cmath>
#include <iostream>

#include "test_main.h"
#include "Matrix.h"

static void testCtorsAndAccess() {
  float a[] = {1, 2, 3, 4, 5, 6};
  Matrix<float> m(a, 2, 3);
  CHECK(m.rows() == 2 && m.cols() == 3);
  CHECK(m.get(0, 0) == 1 && m.get(0, 2) == 3 && m.get(1, 0) == 4 && m.get(1, 2) == 6);

  Matrix<float> c(m);            // copy ctor
  CHECK(c.get(1, 1) == 5);
  Matrix<float> mv = m;          // copy assign
  CHECK(mv.get(1, 2) == 6);
  Matrix<float> moved = std::move(mv);
  CHECK(moved.get(1, 2) == 6 && mv.rows() == 0);
}

static void testRowsCols() {
  float a[] = {1, 2, 3, 4, 5, 6};
  Matrix<float> m(a, 2, 3);
  Vec<float> r0 = m.getrow(0);
  CHECK(r0.len() == 3 && r0[0] == 1 && r0[2] == 3);
  Vec<float> c1 = m.getcol(1);
  CHECK(c1.len() == 2 && c1[0] == 2 && c1[1] == 5);
  Vec<float> nr(3); nr[0] = 7; nr[1] = 8; nr[2] = 9;
  m.setrow(1, nr);
  CHECK(m.get(1, 0) == 7 && m.get(1, 2) == 9);
  Vec<float> nc(2); nc[0] = 9; nc[1] = 8;
  m.setcol(0, nc);
  CHECK(m.get(0, 0) == 9 && m.get(1, 0) == 8);
  // incompatible row/col sizes must throw
  Vec<float> wrongRow(2);
  CHECK_THROWS(m.setrow(0, wrongRow));
  Vec<float> wrongCol(3);
  CHECK_THROWS(m.setcol(0, wrongCol));
}

static void testScalarOps() {
  float a[] = {1, 2, 3, 4};
  Matrix<float> m(a, 2, 2);
  Matrix<float> s = m * 2.0f;
  CHECK(s.get(0, 0) == 2 && s.get(0, 1) == 4 && s.get(1, 0) == 6 && s.get(1, 1) == 8);
  Matrix<float> p = m + 1.0f;
  Matrix<float> q = s / 2.0f;
  CHECK(q.get(0, 0) == 1 && q.get(0, 1) == 2 && q.get(1, 0) == 3 && q.get(1, 1) == 4);
  Matrix<float> mi = m - 1.0f;
  CHECK(mi.get(0, 0) == 0);
  m += 2.0f; m -= 2.0f; m *= 3.0f; m /= 3.0f;
  CHECK(m.get(1, 1) == 4);
  CHECK_THROWS(m / 0.0f);
}

static void testElementwise() {
  float a[] = {1, 2, 3, 4};
  float b[] = {2, 2, 4, 4};
  Matrix<float> m(a, 2, 2), n(b, 2, 2);
  Matrix<float> s = m + n;
  CHECK(s.get(0, 0) == 3 && s.get(1, 1) == 8);
  Matrix<float> d = n - m;
  CHECK(d.get(0, 0) == 1 && d.get(1, 1) == 0);
  Matrix<float> pr = m * n;
  CHECK(pr.get(0, 0) == 2 && pr.get(1, 1) == 16);
  Matrix<float> q = n / m;
  CHECK_NEAR(q.get(1, 0), 4.0 / 3.0, 1e-6);
  // mismatched dims must throw instead of reading OOB
  float c[] = {1, 2, 3};
  Matrix<float> bad(c, 3, 1);
  CHECK_THROWS(m + bad);
  CHECK_THROWS(m += bad);
  CHECK_THROWS(m - bad);
  CHECK_THROWS(m -= bad);
  CHECK_THROWS(m * bad);
  CHECK_THROWS(m *= bad);
  CHECK_THROWS(m / bad);
  CHECK_THROWS(m /= bad);
}

static void testDotTransposeNormalize() {
  float a[] = {1, 2, 3, 4, 5, 6};
  Matrix<float> m(a, 2, 3);
  Vec<float> v(3); v[0] = 1; v[1] = 1; v[2] = 1;
  Vec<float> mv = m.dot(v);
  CHECK(mv.len() == 2 && mv[0] == 6 && mv[1] == 15);

  Matrix<float> t = m.transpose();
  CHECK(t.rows() == 3 && t.cols() == 2 && t.get(0, 0) == 1 && t.get(2, 1) == 6);

  Matrix<float> i2(2, 3);
  i2.set(0, 0, 1); i2.set(1, 0, 2); i2.set(2, 0, 3); // column
  Matrix<float> A = m;
  A.set(0, 0, 0); A.set(0, 1, 0); A.set(0, 2, 0); // zero column
  Matrix<float> nz = A.normalize();
  // zero column stays zero instead of NaN (regression: div by max-min == 0)
  for (u32 c = 0; c < 3; c++)
    for (u32 r = 0; r < 2; r++) {
      double e = nz.get(r, c);
      CHECK(e == e); // not NaN
      CHECK(e >= 0.0 && e <= 1.0);
    }
}

static void testEye() {
  Matrix<float> m(3, 3);
  m.eye(3);
  for (u32 i = 0; i < 3; i++)
    for (u32 j = 0; j < 3; j++)
      CHECK(m.get(i, j) == ((i == j) ? 1 : 0));
}

static void testLinearSolvers() {
  float a[] = {4, 2, 1, 3};  // invertible 2x2
  Matrix<float> m(a, 2, 2);
  Vec<float> b(2); b[0] = 5; b[1] = 6;
  // 4x + 2y = 5, 1x + 3y = 6  =>  x = 3/10, y = 19/10
  Vec<float> x1 = m.solve_1(b);
  CHECK_NEAR(x1[0], 0.3, 1e-4);
  CHECK_NEAR(x1[1], 1.9, 1e-4);

  Vec<float> x2 = m.solve_2(b);
  CHECK_NEAR(x2[0], 0.3, 1e-4);
  CHECK_NEAR(x2[1], 1.9, 1e-4);

  Matrix<float> inv1 = m.inverse_1();
  CHECK_NEAR(inv1.get(0, 0) * m.get(0, 0) + inv1.get(0, 1) * m.get(1, 0), 1.0, 1e-3);
  Matrix<float> inv2 = m.inverse_2();
  CHECK_NEAR(inv2.get(0, 0) * m.get(0, 0) + inv2.get(0, 1) * m.get(1, 0), 1.0, 1e-4);

  CHECK_NEAR(m.determinant_1(), 4.0 * 3 - 2.0 * 1, 1e-3);
  CHECK_NEAR(m.determinant_2(), 4.0 * 3 - 2.0 * 1, 1e-4);
}

static void testQREigen() {
  // QR iteration must recover a diagonal matrix's eigenvalues
  float dia[] = {4, 0, 0, 1};
  Matrix<float> m(dia, 2, 2);
  Vec<double> eig;
  Matrix<double> vec;
  m.eigenDecomposition(eig, vec, 100, 1e-10);
  CHECK_NEAR(eig[0], 4.0, 1e-5);
  CHECK_NEAR(eig[1], 1.0, 1e-5);
}

int main() {
  testCtorsAndAccess();
  testRowsCols();
  testScalarOps();
  testElementwise();
  testDotTransposeNormalize();
  testEye();
  testLinearSolvers();
  testQREigen();
  return test_report();
}
