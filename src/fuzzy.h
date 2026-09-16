#ifndef __FUZZY_H__
#define __FUZZY_H__

#include <cmath>
#include <vector>
#include <limits>
#include <random>

// Euclidean distance between a point and a cluster center
template <typename S, typename T>
double distance(const Matrix<S>& pts, u32 i, const Vec<T>& center) {
    double dist2 = 0.0;
    for (u32 d = 0; d < pts.cols(); ++d) {
        double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(center[d]);
        dist2 += diff * diff;
    }
    return std::sqrt(dist2);
}

template <typename S, typename T>
void fuzzy(const Matrix<S> &pts, const int K, Vec<u8> &labels, Matrix<T> &clusters) {
    const u32 n   = pts.rows();
    const u32 dim = pts.cols();
    const double m = 2.0;   // fuzziness parameter (>1)
    const double tol = 1e-4;
    const int max_iter = 100;

    // Membership matrix U (n x K)
    Matrix<double> U(n, K);

    // Random initialization of memberships
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (u32 i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = 0; k < K; ++k) {
            double val = dist(gen);
            U.set(i, k, val);
            sum += val;
        }
        // normalize
        for (int k = 0; k < K; ++k) {
            U.set(i, k, U.get(i, k) / sum);
        }
    }

    Vec<double> prev_centers(K * dim);
    clusters.init(K, dim, 0);

    for (int iter = 0; iter < max_iter; ++iter) {
        // Step 1: update cluster centers
        for (int k = 0; k < K; ++k)
            for (u32 d = 0; d < dim; ++d)
                prev_centers[k * dim + d] = clusters.get(k, d);
        for (int k = 0; k < K; ++k) {
            for (u32 d = 0; d < dim; ++d) {
                double num = 0.0, den = 0.0;
                for (u32 i = 0; i < n; ++i) {
                    double uik = std::pow(U.get(i, k), m);
                    num += uik * pts.get(i, d);
                    den += uik;
                }
                clusters.set(k, d, (den > 0.0) ? num / den : 0.0);
            }
        }

        // Step 2: update memberships
        for (u32 i = 0; i < n; ++i) {
            for (int k = 0; k < K; ++k) {
                double dist_ik = 0.0;
                for (u32 d = 0; d < dim; ++d) {
                    double diff = pts.get(i, d) - clusters.get(k, d);
                    dist_ik += diff * diff;
                }
                dist_ik = std::sqrt(dist_ik);
                if (dist_ik < 1e-6) {
                    // point coincides with this center: give it full membership
                    for (int kk = 0; kk < K; ++kk) U.set(i, kk, (kk == k) ? 1.0 : 0.0);
                    continue;
                }
                double denom = 0.0;
                for (int j = 0; j < K; ++j) {
                    double dist_ij = 0.0;
                    for (u32 d = 0; d < dim; ++d) {
                        double diff = pts.get(i, d) - clusters.get(j, d);
                        dist_ij += diff * diff;
                    }
                    dist_ij = std::sqrt(dist_ij);
                    denom += std::pow(dist_ik / (dist_ij + 1e-6), 2.0 / (m - 1.0));
                }
                U.set(i, k, 1.0 / denom);
            }
        }

        // Convergence check: max center movement this iteration
        double max_change = 0.0;
        for (int k = 0; k < K; ++k) {
            for (u32 d = 0; d < dim; ++d) {
                double diff = std::abs(clusters.get(k, d) - prev_centers[k * dim + d]);
                max_change = std::max(max_change, diff);
            }
        }
        if (iter > 0 && max_change < tol) break;
    }

    // Assign hard labels based on max membership
    labels.init(n, 0);
    for (u32 i = 0; i < n; ++i) {
        int best_k = 0;
        double best_u = U.get(i, 0);
        for (int k = 1; k < K; ++k) {
            double uik = U.get(i, k);
            if (uik > best_u) {
                best_u = uik;
                best_k = k;
            }
        }
        labels[i] = static_cast<u8>(best_k);
    }
}

#endif // __FUZZY_H__
