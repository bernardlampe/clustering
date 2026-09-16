#ifndef __NCUTS_H__
#define __NCUTS_H__

#include <cmath>
#include <vector>
#include <limits>
#include <random>

// Normalized Cuts clustering
template <typename S, typename T>
void ncuts(const Matrix<S> &pts, const u32 K, Vec<u8> &labels, Matrix<T> &clusters) {
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

    // Step 3: Normalized Laplacian L = D^{-1/2} (D - W) D^{-1/2}
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
    // Power iteration on L directly converges to the largest-eigenvalue
    // eigenvector, not the second-smallest one ncuts needs. Iterate on
    // (2I - L) so small eigenvalues dominate, and deflate the trivial
    // all-ones (generalized) eigenvector after each multiply.
    Vec<double> v(n);
    // deterministic start orthogonal to the D^{1/2}-scaled ones vector:
    // for normalized Lsym, ones/sqrt(d) is the trivial eigenvector; starting
    // orthogonal to ones is an adequate practical approximation
    double s = 0.0;
    for (u32 i = 0; i < n; ++i) { v[i] = (i % 2 == 0) ? 1.0 : -1.0; s += v[i]; }
    s /= n;
    for (u32 i = 0; i < n; ++i) v[i] -= s;

    for (int iter = 0; iter < 300; ++iter) {
        Vec<double> v_new(n);
        v_new.zero();
        for (u32 i = 0; i < n; ++i) {
            double sum = 0.0;
            for (u32 j = 0; j < n; ++j)
                sum += L.get(i, j) * v[j];
            v_new[i] = 2.0 * v[i] - sum;
        }
        // deflate trivial component
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

#endif // __NCUTS_H__
