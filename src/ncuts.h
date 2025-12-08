#ifndef __NCUTS_H__
#define __NCUTS_H__

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

// Normalized Cuts clustering
template <typename S, typename T>
void ncuts(const Matrix<S> &pts, const u32 K, Vec<u8> &labels, Matrix<T> &clusters) {
    const u32 n = pts.rows();
    const u32 dim = pts.cols();

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

    // Step 3: Normalized Laplacian L = D^{-1/2} (D - W) D^{-1/2}
    Matrix<double> L(n, n);
    for (u32 i = 0; i < n; ++i) {
        for (u32 j = 0; j < n; ++j) {
            if (i == j && D.get(i, i) > 0) {
                L.set(i, j, 1.0);
            } else if (D.get(i, i) > 0 && D.get(j, j) > 0) {
                L.set(i, j, -W.get(i, j) / std::sqrt(D.get(i, i) * D.get(j, j)));
            }
        }
    }

    // Step 4: Eigen decomposition (simplified: use power iteration for 2nd smallest eigenvector)
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

#endif // __NCUTS_H__
