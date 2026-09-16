#ifndef __SPECTRAL_H__
#define __SPECTRAL_H__

#include <cmath>
#include <vector>
#include <limits>
#include <random>

// Spectral clustering
template <typename S, typename T>
void spectral(const Matrix<S> &pts, const u32 K, Vec<u8> &labels, Matrix<T> &clusters) {
    const u32 n = pts.rows();

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

    // Step 4: Fiedler vector via shifted power iteration.
    // Plain power iteration on L converges to the LARGEST eigenvalue's
    // eigenvector, which is useless for the embedding; we need the smallest
    // non-trivial one. Iterate on (2I - L) so that eigenvalue lambda -> 2 -
    // lambda: the smallest lambda becomes the largest 2-lambda. After each
    // multiply, deflate the constant vector (the trivial eigenpair of L),
    // otherwise the iteration converges back to it.
    Vec<double> v(n);
    // deterministic random-ish start, orthogonal to the all-ones vector
    double s = 0.0;
    for (u32 i = 0; i < n; ++i) { v[i] = (i % 2 == 0) ? 1.0 : -1.0; s += v[i]; }
    s /= n;
    for (u32 i = 0; i < n; ++i) v[i] -= s; // start orthogonal to ones

    for (int iter = 0; iter < 300; ++iter) {
        Vec<double> v_new(n);
        v_new.zero();
        for (u32 i = 0; i < n; ++i) {
            double sum = 0.0;
            for (u32 j = 0; j < n; ++j)
                sum += L.get(i, j) * v[j];
            v_new[i] = 2.0 * v[i] - sum;
        }
        // deflate trivial eigenvector (all-ones) component
        double mean = 0.0;
        for (u32 i = 0; i < n; ++i) mean += v_new[i];
        mean /= n;
        for (u32 i = 0; i < n; ++i) v_new[i] -= mean;

        // normalize
        double norm = 0.0;
        for (u32 i = 0; i < n; ++i) norm += v_new[i] * v_new[i];
        norm = std::sqrt(norm);
        if (norm < 1e-12) break;
        for (u32 i = 0; i < n; ++i) v[i] = v_new[i] / norm;
    }

    // Step 5: Use Fiedler vector as embedding
    Matrix<T> embed(n, 1);
    for (u32 i = 0; i < n; ++i) embed.set(i, 0, static_cast<T>(v[i]));

    // Step 6: Cluster embedding with k-means
    kmeans(embed, K, labels, clusters);
}

#endif // __SPECTRAL_H__
