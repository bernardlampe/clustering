#ifndef __DBSCAN_H__
#define __DBSCAN_H__

#include <cmath>
#include <vector>

// Helper: Euclidean distance
template <typename S>
double distance(const Matrix<S>& pts, u32 i, u32 j) {
    double dist2 = 0.0;
    for (u32 d = 0; d < pts.cols(); ++d) {
        double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(pts.get(j, d));
        dist2 += diff * diff;
    }
    return std::sqrt(dist2);
}

// Find neighbors of point i within radius
template <typename S>
std::vector<u32> regionQuery(const Matrix<S>& pts, u32 i, double radius) {
    std::vector<u32> neighbors;
    for (u32 j = 0; j < pts.rows(); ++j) {
        if (distance(pts, i, j) <= radius) {
            neighbors.push_back(j);
        }
    }
    return neighbors;
}

template <typename S, typename T>
void dbscan(const Matrix<S> &pts, const float radius, const int minpts,
            Vec<u8> &labels, Matrix<T> &clusters) {
    const u32 n   = pts.rows();
    const u32 dim = pts.cols();

    // Special label values
    const u8 UNVISITED = 255; // internal marker
    const u8 NOISE     = 254; // internal marker

    labels.init(n, UNVISITED);
    int cluster_id = 0;

    for (u32 i = 0; i < n; ++i) {
        if (labels[i] != UNVISITED) continue; // already processed

        auto neighbors = regionQuery(pts, i, radius);

        if (neighbors.size() < static_cast<u32>(minpts)) {
            labels[i] = NOISE;
            continue;
        }

        // Start new cluster
        labels[i] = static_cast<u8>(cluster_id);

        std::vector<u32> seeds = neighbors;
        for (size_t idx = 0; idx < seeds.size(); ++idx) {
            u32 p = seeds[idx];

            if (labels[p] == NOISE) {
                labels[p] = static_cast<u8>(cluster_id); // noise becomes border
            }
            if (labels[p] != UNVISITED) continue;

            labels[p] = static_cast<u8>(cluster_id);
            auto p_neighbors = regionQuery(pts, p, radius);
            if (p_neighbors.size() >= static_cast<u32>(minpts)) {
                seeds.insert(seeds.end(), p_neighbors.begin(), p_neighbors.end());
            }
        }
        cluster_id++;
    }

    // Compute cluster centers (mean of points in each cluster)
    clusters.init(cluster_id, dim, 0);
    std::vector<int> counts(cluster_id, 0);

    for (u32 i = 0; i < n; ++i) {
        int cid = labels[i];
        if (cid >= 0 && cid < cluster_id) {
            counts[cid]++;
            for (u32 d = 0; d < dim; ++d) {
                clusters.set(cid, d, clusters.get(cid, d) + pts.get(i, d));
            }
        }
    }

    for (int c = 0; c < cluster_id; ++c) {
        if (counts[c] > 0) {
            for (u32 d = 0; d < dim; ++d) {
                clusters.set(c, d, clusters.get(c, d) / counts[c]);
            }
        }
    }
}


#endif // __DBSCAN_H__
