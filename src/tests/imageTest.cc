#include <cstdio>
#include <cmath>
#include <string>

#include "test_main.h"
#include "Image.h"

static void testReadWriteRoundtrip() {
  // u8 pgm
  {
    Image<u8> im(4, 7);
    for (u32 i = 0; i < 4 * 7; i++) im[i] = (u8)(i * 3);
    im.writeToFile("/tmp/t_pgm.pgm");

    Image<u8> back;
    back.readFromFile("/tmp/t_pgm.pgm");
    CHECK(back.height() == 4 && back.width() == 7);
    bool ok = true;
    for (u32 i = 0; i < 4 * 7; i++) ok = ok && back[i] == (u8)(i * 3);
    CHECK(ok); // raster not shifted by header newline
  }
  // RGB ppm
  {
    Image<RGB_t> im(3, 5);
    for (u32 i = 0; i < 3 * 5; i++) {
      im[i][0] = (u8)(i % 256);
      im[i][1] = (u8)(255 - (i % 256));
      im[i][2] = (u8)(128);
    }
    im.writeToFile("/tmp/t_ppm.ppm");

    Image<RGB_t> back;
    back.readFromFile("/tmp/t_ppm.ppm");
    CHECK(back.height() == 3 && back.width() == 5);
    bool ok = true;
    for (u32 i = 0; i < 3 * 5; i++) {
      ok = ok && back[i][0] == im[i][0] && back[i][1] == im[i][1] && back[i][2] == im[i][2];
    }
    CHECK(ok);
  }
}

static void testReaderRejectsWrongMagic() {
  {
    FILE *f = fopen("/tmp/t_bad.pgm", "wb");
    fputs("P4\n2 2\n255\n", f);
    for (int i = 0; i < 4; i++) fputc(0, f);
    fclose(f);
    Image<u8> im;
    CHECK_THROWS(im.readFromFile("/tmp/t_bad.pgm"));
  }
  {
    FILE *f = fopen("/tmp/t_bad.ppm", "wb");
    fputs("P5\n2 2\n255\n", f);
    for (int i = 0; i < 4; i++) fputc(0, f);
    fclose(f);
    Image<RGB_t> im;
    CHECK_THROWS(im.readFromFile("/tmp/t_bad.ppm"));
  }
}

static void testConvolveSeparableMatchesReference() {
  // regression: 1-D convolve previously wrote transposed pixels; for a
  // non-square image results were scatter-corrupted
  const u32 H = 2, W = 3;
  Image<float> im(H, W);
  float v = 0;
  for (u32 h = 0; h < H; h++)
    for (u32 w = 0; w < W; w++)
      im.set(h, w, ++v); // {1 2 3; 4 5 6}

  float k[3] = {0.25f, 0.5f, 0.25f};
  // reference computed independently (python) for this kernel and image
  const float expected[H * W] = {1.3125f, 2.25f, 2.0625f, 1.875f, 3.0f, 2.625f};
  im.convolve(k, 3);
  bool ok = true;
  for (u32 h = 0; h < H; h++)
    for (u32 w = 0; w < W; w++)
      ok = ok && std::fabs(im.get(h, w) - expected[h * W + w]) < 1e-4;
  CHECK(ok);
}

static void testConvolve2DIdentityAndBox() {
  const u32 H = 4, W = 5;
  Image<float> im(H, W);
  for (u32 i = 0; i < H * W; i++) im[i] = (float)(i % 11) + 1.0f;

  // identity kernel must not change pixels
  float id[9] = {0, 0, 0, 0, 1, 0, 0, 0, 0};
  Image<float> before = im;
  im.convolve(id, 3, 3);
  bool same = true;
  for (u32 i = 0; i < H * W; i++) same = same && im[i] == before[i];
  CHECK(same);

  // 3x3 box kernel: interior pixels average the full 9-window, border
  // pixels use a truncated window (out-of-bounds neighbors are skipped,
  // not zero-padded): interior 9, edges 6, corners 4
  float box[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  Image<float> b(H, W);
  for (u32 i = 0; i < H * W; i++) b[i] = 1.0f;
  b.convolve(box, 3, 3);
  CHECK(b.get(0, 0) == 4.0f);          // corner
  CHECK(b.get(1, 0) == 6.0f);          // left edge
  CHECK(b.get(1, 1) == 9.0f);          // interior
  CHECK(b.get(3, 4) == 4.0f);          // far corner
}

static void testCreatepts() {
  Image<u8> g(2, 3);
  for (u32 i = 0; i < 6; i++) g[i] = (u8)(i + 10);
  Matrix<float> pts = g.createpts();
  CHECK(pts.rows() == 6 && pts.cols() == 1);
  CHECK(pts.get(0, 0) == 10 && pts.get(5, 0) == 15);

  Image<RGB_t> c(2, 2);
  for (u32 i = 0; i < 4; i++) {
    c[i][0] = (u8)(i * 1);
    c[i][1] = (u8)(i * 2);
    c[i][2] = (u8)(i * 3);
  }
  Matrix<float> cpts = c.createpts();
  CHECK(cpts.rows() == 4 && cpts.cols() == 3);
  CHECK(cpts.get(3, 0) == 3 && cpts.get(3, 1) == 6 && cpts.get(3, 2) == 9);
}

static void testArithmeticOperators() {
  Image<u8> a(2, 2), b(2, 2);
  for (u32 i = 0; i < 4; i++) { a[i] = 10; b[i] = 4; }
  Image<u8> s = a + b;
  CHECK(s[0] == 14 && s[3] == 14);
  Image<u8> d = a - b;
  CHECK(d[0] == 6);
  Image<u8> p = a * b;
  CHECK(p[0] == 40);
}

int main() {
  testReadWriteRoundtrip();
  testReaderRejectsWrongMagic();
  testConvolveSeparableMatchesReference();
  testConvolve2DIdentityAndBox();
  testCreatepts();
  testArithmeticOperators();
  return test_report();
}
