#ifndef __KMEANS_H__
#define __KMEANS_H__

#include<limits>

template <typename S, typename T>
void kmeans(const Matrix<S> &pts, const u32 nclusters, Vec<u8> &labels, Matrix<T> &clusters) {
  bool change;
  double mind, d;
  u8 minl;
  Vec<S> pt;
  Vec<float> mean;
  u32 dim = pts.cols();
  u32 npts = pts.rows();
  u32 count;

  // check nclusters
  if (nclusters <= 0)
    throw Exception("kmeans must have nclusters > 0");

  // assign each pts to a random label  
  labels.init(npts);
  for(u32 i = 0; i < labels.len(); i++) labels[i] = i % nclusters;

  // kmeans algo loop
  mean.init(dim);
  clusters.init(nclusters, dim);
  do {
    // assume perfect labeling
    change = false;

    // compute means of each cluster
    for(u32 i = 0; i < nclusters; i++) {
        mean.zero();
        count = 0;

        for(u32 j = 0; j < npts; j++) {
            if (labels[j] == i) {
                mean += pts.getrow(j); // add point to mean
                count++;
            }
        }

        if (count == 0) {
            // empty cluster: re-seed it with a random point instead of
            // dividing by zero, which poisons the centers with NaN
            u32 seed = rand() % npts;
            clusters.setrow(i, pts.getrow(seed));
            continue;
        }

        // put mean in cluster matrix
        mean /= (double)count;
        clusters.setrow(i, mean);
    }

    // update labels for each pt to nearest cluster, flag if any change
    for(u32 i = 0; i < npts; i++) {
        pt = pts.getrow(i);
        mind = std::numeric_limits<double>::max();
        for(u32 j = 0; j < nclusters; j++) {
            d = (pt - clusters.getrow(j)).norm();
            if (d < mind) {
                mind = d;
                minl = j;
            }
        }

        if (minl != labels[i]) {
            change = true;
            labels[i] = minl;
        }
    }
  } while (change);
}

#endif // __KMEANS_H__
