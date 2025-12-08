#ifndef __SPECTRAL_H__
#define __SPECTRAL_H__

#include <cmath>
#include <vector>
#include <limits>
#include <random>

// Gaussian similarity function
template <typename S>
double similarity(const Matrix<S>& pts, u32 i, u32 j, double sigma = 1.0) {
    double dist2 = 0.0;
    for (u32 d = 0; d < pts.cols(); ++d) {
        double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(pts.get(j, d));
        dist2 += diff * diff;
    }
    return std::exp(-dist2 / (2 * sigma * sigma));
}

// Simple k-means clustering
template <typename T>
void kmeans(const Matrix<T>& embed, int K, Vec<u8>& labels, Matrix<T>& clusters) {
    const u32 n = embed.rows();
    const u32 dim = embed.cols();

    clusters.init(K, dim);
    labels.init(n, 0);

    // Random initialization of cluster centers
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<u32> dist(0, n - 1);
    for (int k = 0; k < K; ++k) {
        u32 idx = dist(gen);
        for (u32 d = 0; d < dim; ++d)
            clusters.set(k, d, embed.get(idx, d));
    }

    bool changed = true;
    int max_iter = 100;
    while (changed && max_iter--) {
        changed = false;

        // Assign labels
        for (u32 i = 0; i < n; ++i) {
            double bestDist = std::numeric_limits<double>::infinity();
            int bestK = 0;
            for (int k = 0; k < K; ++k) {
                double dist2 = 0.0;
                for (u32 d = 0; d < dim; ++d) {
                    double diff = embed.get(i, d) - clusters.get(k, d);
                    dist2 += diff * diff;
                }
                if (dist2 < bestDist) {
                    bestDist = dist2;
                    bestK = k;
                }
            }
            if (labels[i] != bestK) {
                labels[i] = bestK;
                changed = true;
            }
        }

        // Update cluster centers
        clusters.init(K, dim, 0);
        std::vector<int> counts(K, 0);
        for (u32 i = 0; i < n; ++i) {
            int cid = labels[i];
            counts[cid]++;
            for (u32 d = 0; d < dim; ++d)
                clusters.set(cid, d, clusters.get(cid, d) + embed.get(i, d));
        }
        for (int k = 0; k < K; ++k) {
            if (counts[k] > 0) {
                for (u32 d = 0; d < dim; ++d)
                    clusters.set(k, d, clusters.get(k, d) / counts[k]);
            }
        }
    }
}

// Spectral clustering
template <typename S, typename T>
void spectral(const Matrix<S> &pts, Vec<u8> &labels, Matrix<T> &clusters) {
    const u32 n = pts.rows();
    const u32 dim = pts.cols();
    const int K = 2; // number of clusters (can be parameterized)

    // Step 1: Build similarity matrix W
    Matrix<double> W(n, n);
    for (u32 i = 0; i < n; ++i) {
        for (u32 j = 0; j < n; ++j) {
            if (i == j) W.set(i, j, 0.0);
            else W.set(i, j, similarity(pts, i, j, 1.0));
        }
    }

    // Step 2: Degree matrix D
    Matrix<double> D(n, n);
    for (u32 i = 0; i < n; ++i) {
        double sum = 0.0;
        for (u32 j = 0; j < n; ++j) sum += W.get(i, j);
        D.set(i, i, sum);
    }

    // Step 3: Normalized Laplacian L = I - D^{-1/2} W D^{-1/2}
    Matrix<double> L(n, n);
    for (u32 i = 0; i < n; ++i) {
        for (u32 j = 0; j < n; ++j) {
            if (i == j) {
                L.set(i, j, 1.0);
            } else if (D.get(i, i) > 0 && D.get(j, j) > 0) {
                L.set(i, j, -W.get(i, j) / std::sqrt(D.get(i, i) * D.get(j, j)));
            }
        }
    }

    // Step 4: Approximate eigenvector (power iteration)
    Vec<double> v(n);
    v.one(); // initial vector
    for (int iter = 0; iter < 100; ++iter) {
        Vec<double> v_new(n);
        for (u32 i = 0; i < n; ++i) {
            double sum = 0.0;
            for (u32 j = 0; j < n; ++j)
                sum += L.get(i, j) * v[j];
            v_new[i] = sum;
        }
        // normalize
        double norm = 0.0;
        for (u32 i = 0; i < n; ++i) norm += v_new[i] * v_new[i];
        norm = std::sqrt(norm);
        for (u32 i = 0; i < n; ++i) v[i] = v_new[i] / norm;
    }

    // Step 5: Use eigenvector as embedding
    Matrix<T> embed(n, 1);
    for (u32 i = 0; i < n; ++i) embed.set(i, 0, static_cast<T>(v[i]));

    // Step 6: Cluster embedding with k-means
    kmeans(embed, K, labels, clusters);
}


#endif // __SPECTRAL_H__
