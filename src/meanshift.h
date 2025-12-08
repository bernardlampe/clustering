#ifndef __MEANSHIFT_H__
#define __MEANSHIFT_H__


#include <cmath>
#include <vector>
#include <limits>

// Gaussian kernel function
template <typename S>
double gaussian_kernel(const Matrix<S> &pts, u32 i, const Vec<S> &center, double bandwidth) {
    double dist2 = 0.0;
    for (u32 d = 0; d < pts.cols(); d++) {
        double diff = pts.get(i, d) - center[d];
        dist2 += diff * diff;
    }
    return std::exp(-dist2 / (2 * bandwidth * bandwidth));
}

template <typename S, typename T>
void meanshift(const Matrix<S> &pts, const double bandwidth, const u32 max_iters, const double tol, Vec<u8> &labels, Matrix<T> &clusters) {
    const u32 n = pts.rows();
    const u32 dim = pts.cols();

    // Initialize shifted points as input points
    Matrix<T> shifted(n, dim);
    for (u32 i = 0; i < n; i++)
        for (u32 d = 0; d < dim; d++)
            shifted.set(i, d, pts.get(i, d));

    // Iteratively shift each point
    for (u32 iter = 0; iter < max_iters; iter++) {
        bool converged = true;
        for (u32 i = 0; i < n; i++) {
            Vec<T> new_center(dim);
            double weight_sum = 0.0;

            // Compute weighted mean of neighbors
            for (u32 j = 0; j < n; j++) {
                double w = gaussian_kernel(pts, j, shifted.getrow(i), bandwidth);
                weight_sum += w;
                for (u32 d = 0; d < dim; d++)
                    new_center[d] += w * pts.get(j, d);
            }

            for (u32 d = 0; d < dim; d++)
                new_center[d] /= weight_sum;

            // Check convergence
            double diff_norm = 0.0;
            for (u32 d = 0; d < dim; d++) {
                double diff = shifted.get(i, d) - new_center[d];
                diff_norm += diff * diff;
                shifted.set(i, d, new_center[d]);
            }
            if (std::sqrt(diff_norm) > tol) converged = false;
        }
        if (converged) break;
    }

    // Identify unique clusters (merge close centers)
    std::vector<Vec<T>> centers;
    labels.init(n, 0);

    for (u32 i = 0; i < n; i++) {
        bool found = false;
        for (u32 c = 0; c < centers.size(); c++) {
            double dist2 = 0.0;
            for (u32 d = 0; d < dim; d++) {
                double diff = shifted.get(i, d) - centers[c][d];
                dist2 += diff * diff;
            }
            if (std::sqrt(dist2) < bandwidth / 2) {
                labels[i] = c;
                found = true;
                break;
            }
        }
        if (!found) {
            Vec<T> new_center(dim);
            for (u32 d = 0; d < dim; d++)
                new_center[d] = shifted.get(i, d);
            centers.push_back(new_center);
            labels[i] = centers.size() - 1;
        }
    }

    // Save clusters into output matrix
    clusters.init(centers.size(), dim);
    for (u32 c = 0; c < centers.size(); c++)
        for (u32 d = 0; d < dim; d++)
            clusters.set(c, d, centers[c][d]);
}

#endif // __MEANSHIFT_H__
