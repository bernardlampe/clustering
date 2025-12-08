#ifndef __EM_H__
#define __EM_H__

#include <cmath>
#include <random>
#include <limits>

// Helper: squared Euclidean distance between a point and a mean
template <typename S, typename T>
double sqdist_point_mean(const Matrix<S>& pts, u32 i, const Vec<T>& mean) {
    double d2 = 0.0;
    for (u32 d = 0; d < pts.cols(); ++d) {
        double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(mean[d]);
        d2 += diff * diff;
    }
    return d2;
}

// Helper: k-means++-style initialization of means
template <typename S, typename T>
void init_means_kpp(const Matrix<S>& pts, Matrix<T>& means, u32 K) {
    const u32 n = pts.rows();
    const u32 dim = pts.cols();
    means.init(K, dim);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<u32> uni_idx(0, n - 1);

    // pick first center randomly
    u32 first = uni_idx(gen);
    for (u32 d = 0; d < dim; ++d) means.set(0, d, static_cast<T>(pts.get(first, d)));

    // distances to nearest chosen center
    std::vector<double> minD2(n, std::numeric_limits<double>::infinity());

    for (u32 k = 1; k < K; ++k) {
        // update min distances
        for (u32 i = 0; i < n; ++i) {
            // compute distance to latest center k-1
            double d2 = 0.0;
            for (u32 d = 0; d < dim; ++d) {
                double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(means.get(k - 1, d));
                d2 += diff * diff;
            }
            if (d2 < minD2[i]) minD2[i] = d2;
        }

        // sample next center proportional to D^2
        double sumD2 = 0.0;
        for (double v : minD2) sumD2 += v;
        std::uniform_real_distribution<double> uni_prob(0.0, sumD2);
        double r = uni_prob(gen);

        u32 chosen = 0;
        double accum = 0.0;
        for (u32 i = 0; i < n; ++i) {
            accum += minD2[i];
            if (accum >= r) { chosen = i; break; }
        }
        for (u32 d = 0; d < dim; ++d) means.set(k, d, static_cast<T>(pts.get(chosen, d)));
    }
}

// EM for Gaussian Mixture Models (diagonal covariance)
template <typename S, typename T>
void em(const Matrix<S> &pts, Vec<u8> &labels, Matrix<T> &clusters) {
    const u32 n = pts.rows();
    const u32 dim = pts.cols();

    // Choose number of clusters (can be adjusted)
    const u32 K = 3;

    // Parameters: mixing coefficients, means (clusters), diagonal variances
    Vec<double> pi(K);           // mixture weights
    Matrix<T> means;             // K x dim
    Matrix<double> var;          // K x dim (diagonal covariances)

    // Responsibilities: n x K
    Matrix<double> resp(n, K);

    // Initialization
    init_means_kpp(pts, means, K);
    clusters = means; // keep clusters in sync with means
    for (u32 k = 0; k < K; ++k) pi[k] = 1.0 / static_cast<double>(K);

    var.init(K, dim, 1.0); // start with unit variances

    const double eps = 1e-6;      // variance floor
    const double tol = 1e-4;      // log-likelihood tolerance
    const u32 max_iter = 100;

    double prev_ll = -std::numeric_limits<double>::infinity();

    for (u32 iter = 0; iter < max_iter; ++iter) {
        // E-step: compute responsibilities r_{ik} proportional to pi_k * N(x_i | mu_k, diag(var_k))
        for (u32 i = 0; i < n; ++i) {
            // compute unnormalized responsibilities
            double sum_r = 0.0;
            for (u32 k = 0; k < K; ++k) {
                // log probability of x_i under component k
                double logp = 0.0;
                double logdet = 0.0;
                for (u32 d = 0; d < dim; ++d) {
                    double sigma2 = std::max(var.get(k, d), eps);
                    double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(means.get(k, d));
                    logp += -0.5 * (diff * diff) / sigma2;
                    logdet += std::log(sigma2);
                }
                logp += -0.5 * (dim * std::log(2.0 * M_PI) + logdet);
                double val = std::exp(logp) * pi[k];
                resp.set(i, k, val);
                sum_r += val;
            }
            // normalize responsibilities across k
            double inv_sum = (sum_r > 0.0) ? (1.0 / sum_r) : 1.0 / static_cast<double>(K);
            for (u32 k = 0; k < K; ++k) {
                resp.set(i, k, resp.get(i, k) * inv_sum);
            }
        }

        // M-step: update pi, means, variances using responsibilities
        // Nk = sum_i r_{ik}
        Vec<double> Nk(K);
        Nk.zero();
        for (u32 k = 0; k < K; ++k) {
            double sum_r = 0.0;
            for (u32 i = 0; i < n; ++i) sum_r += resp.get(i, k);
            Nk[k] = sum_r;
            pi[k] = sum_r / static_cast<double>(n);
        }

        // Update means
        for (u32 k = 0; k < K; ++k) {
            for (u32 d = 0; d < dim; ++d) {
                double num = 0.0;
                for (u32 i = 0; i < n; ++i) {
                    num += resp.get(i, k) * static_cast<double>(pts.get(i, d));
                }
                double mu = (Nk[k] > 0.0) ? (num / Nk[k]) : static_cast<double>(means.get(k, d));
                means.set(k, d, static_cast<T>(mu));
            }
        }

        // Update variances (diagonal)
        for (u32 k = 0; k < K; ++k) {
            for (u32 d = 0; d < dim; ++d) {
                double num = 0.0;
                for (u32 i = 0; i < n; ++i) {
                    double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(means.get(k, d));
                    num += resp.get(i, k) * diff * diff;
                }
                double sigma2 = (Nk[k] > 0.0) ? (num / Nk[k]) : var.get(k, d);
                var.set(k, d, std::max(sigma2, eps));
            }
        }

        // Compute log-likelihood
        double ll = 0.0;
        for (u32 i = 0; i < n; ++i) {
            double sum_comp = 0.0;
            for (u32 k = 0; k < K; ++k) {
                double logp = 0.0;
                double logdet = 0.0;
                for (u32 d = 0; d < dim; ++d) {
                    double sigma2 = std::max(var.get(k, d), eps);
                    double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(means.get(k, d));
                    logp += -0.5 * (diff * diff) / sigma2;
                    logdet += std::log(sigma2);
                }
                logp += -0.5 * (dim * std::log(2.0 * M_PI) + logdet);
                sum_comp += pi[k] * std::exp(logp);
            }
            ll += std::log(std::max(sum_comp, 1e-300));
        }

        // Check convergence
        if (std::abs(ll - prev_ll) < tol) break;
        prev_ll = ll;
    }

    // Assign labels based on maximum responsibility
    labels.init(n, 0);
    for (u32 i = 0; i < n; ++i) {
        u32 best_k = 0;
        double best_r = resp.get(i, 0);
        for (u32 k = 1; k < K; ++k) {
            double r = resp.get(i, k);
            if (r > best_r) { best_r = r; best_k = k; }
        }
        labels[i] = static_cast<u8>(best_k);
    }

    // Output cluster centers
    clusters.init(K, dim);
    for (u32 k = 0; k < K; ++k)
        for (u32 d = 0; d < dim; ++d)
            clusters.set(k, d, means.get(k, d));
}

#endif // __EM_H__
