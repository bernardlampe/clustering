#include <cmath>
#include <cstdio>
#include <vector>

#include "test_main.h"
#include "cluster.h"
#include "Matrix.h"
#include "Vec.h"
#include "types.h"

// Three well-separated 2-D blobs plus deterministic construction: cluster
// structure is unambiguous, so every algorithm must recover >= 2 clusters and
// label all points in the same blob mostly alike. Points are NOT normalized
// here beyond what cluster() does internally (it min-max normalizes).

static Matrix<float> makeThreeBlobs() {
  const u32 per = 60; // points per blob
  Matrix<float> pts(per * 3, 2);
  // blob centers in [0,1]^2 with margin >> spread
  const float cx[3] = {0.15f, 0.85f, 0.5f};
  const float cy[3] = {0.15f, 0.85f, 0.85f};
  u32 idx = 0;
  for (u32 b = 0; b < 3; b++) {
    for (u32 i = 0; i < per; i++) {
      // golden-ratio low-discrepancy jitter: 60 unique points per blob,
      // spread 0.01 << blob separation 0.3, deterministic across platforms
      float dx = (0.6180339887f * i - std::floor(0.6180339887f * i) - 0.5f) * 0.02f;
      float dy = (0.3819660113f * i - std::floor(0.3819660113f * i) - 0.5f) * 0.02f;
      pts.set(idx, 0, cx[b] + dx);
      pts.set(idx, 1, cy[b] + dy);
      idx++;
    }
  }
  return pts;
}

static bool labelAgreesWithBlob(const Vec<u8> &labels, u32 per) {
  // labels[i] must be 1:1 with blob index up to relabeling; check that the
  // label set of each blob is dominant (>= 50% of the blob shares one label)
  for (u32 b = 0; b < 3; b++) {
    std::vector<int> count(256, 0);
    for (u32 i = b * per; i < (b + 1) * per; i++) count[labels[i]]++;
    int best = 0;
    for (int c = 0; c < 256; c++) best = std::max(best, count[c]);
    if (best < per / 2) return false;
  }
  return true;
}

static void runAlgo(Params_t params, const Matrix<float> &pts, const char *name,
                    u32 expectMinClusters) {
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(params, pts, labels, clusters);
  CHECK(labels.len() == pts.rows());
  CHECK(clusters.rows() >= expectMinClusters);
  // every point gets a label; cluster centers are finite
  for (u32 i = 0; i < labels.len(); i++) {
    if (labels[i] == 254) continue; // dbscan noise
    CHECK(labels[i] < clusters.rows());
  }
  for (u32 c = 0; c < clusters.rows(); c++)
    for (u32 d = 0; d < clusters.cols(); d++) {
      double v = clusters.get(c, d);
      CHECK(v == v); // not NaN
    }
}

static void testKmeans() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "kmeans"; // legacy dispatch key must still work
  p["nclusters"] = "3";
  runAlgo(p, pts, "kmeans", 1);
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(p, pts, labels, clusters);
  CHECK(labelAgreesWithBlob(labels, 60));
}

static void testDbscan() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "dbscan";
  p["radius"] = "0.1";
  p["minpts"] = "5";
  runAlgo(p, pts, "dbscan", 1);
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(p, pts, labels, clusters);
  // dense separated blobs: (mostly) no noise, and blob-consistent labels
  u32 noise = 0;
  for (u32 i = 0; i < labels.len(); i++)
    if (labels[i] == 254 || labels[i] == 255) noise++;
  CHECK(noise == 0);
  CHECK(labelAgreesWithBlob(labels, 60));
}

static void testEm() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "em";
  p["nclusters"] = "3";
  runAlgo(p, pts, "em", 1);
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(p, pts, labels, clusters);
  CHECK(labelAgreesWithBlob(labels, 60));
}

static void testFuzzy() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "fuzzy";
  p["nclusters"] = "3";
  runAlgo(p, pts, "fuzzy", 1);
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(p, pts, labels, clusters);
  CHECK(labelAgreesWithBlob(labels, 60));
}

static void testMeanshift() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "meanshift";
  p["bandwidth"] = "0.15";
  p["max_iters"] = "50";
  p["tol"] = "0.001";
  runAlgo(p, pts, "meanshift", 1);
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(p, pts, labels, clusters);
  CHECK(clusters.rows() >= 2); // close to k clusters found, never fewer than 2
  CHECK(labelAgreesWithBlob(labels, 60));
}

static void testSpectral() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "spectral";
  p["nclusters"] = "3";
  runAlgo(p, pts, "spectral", 1);
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(p, pts, labels, clusters);
  CHECK(labelAgreesWithBlob(labels, 60));
}

static void testNcuts() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "ncuts";
  p["nclusters"] = "3";
  runAlgo(p, pts, "ncuts", 1);
  Vec<u8> labels;
  Matrix<float> clusters;
  cluster(p, pts, labels, clusters);
  CHECK(labelAgreesWithBlob(labels, 60));
}

static void testUnknownAlgoThrows() {
  Matrix<float> pts = makeThreeBlobs();
  Params_t p;
  p["c"] = "nosuchalgo";
  Vec<u8> labels;
  Matrix<float> clusters;
  CHECK_THROWS(cluster(p, pts, labels, clusters));
}

static void testMissingParamsThrow() {
  Matrix<float> pts = makeThreeBlobs();
  Vec<u8> labels;
  Matrix<float> clusters;
  {
    Params_t p;
    p["c"] = "kmeans";
    CHECK_THROWS(cluster(p, pts, labels, clusters)); // no nclusters
  }
  {
    Params_t p;
    p["c"] = "dbscan";
    p["minpts"] = "5";
    CHECK_THROWS(cluster(p, pts, labels, clusters)); // no radius
  }
}

int main() {
  testKmeans();
  testDbscan();
  testEm();
  testFuzzy();
  testMeanshift();
  testSpectral();
  testNcuts();
  testUnknownAlgoThrows();
  testMissingParamsThrow();
  return test_report();
}
