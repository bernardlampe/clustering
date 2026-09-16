A C++ implementation of unsupervised clustering algorithms applied to image segmentation.

## Build

    make -C src            # build `src/cluster`
    make -C src/tests test # build and run the test suite

## Results

| Dataset | kmeans | dbscan | em | fuzzy | meanshift | spectral | ncuts |
|---------|--------|--------|----|-------|-----------|----------|-------|
| cluster0 | ![`kmeans_0`](images/kmeans_0.png) | ![`dbscan_0`](images/dbscan_0.png) | ![`em_0`](images/em_0.png) | ![`fuzzy_0`](images/fuzzy_0.png) | ![`meanshift_0`](images/meanshift_0.png) | ![`spectral_0`](images/spectral_0.png) | ![`ncuts_0`](images/ncuts_0.png) |
| cluster1 | ![`kmeans_1`](images/kmeans_1.png) | ![`dbscan_1`](images/dbscan_1.png) | ![`em_1`](images/em_1.png) | ![`fuzzy_1`](images/fuzzy_1.png) | ![`meanshift_1`](images/meanshift_1.png) | ![`spectral_1`](images/spectral_1.png) | ![`ncuts_1`](images/ncuts_1.png) |
| cluster2 | ![`kmeans_2`](images/kmeans_2.png) | ![`dbscan_2`](images/dbscan_2.png) | ![`em_2`](images/em_2.png) | ![`fuzzy_2`](images/fuzzy_2.png) | ![`meanshift_2`](images/meanshift_2.png) | ![`spectral_2`](images/spectral_2.png) | ![`ncuts_2`](images/ncuts_2.png) |
| cluster3 | ![`kmeans_3`](images/kmeans_3.png) | ![`dbscan_3`](images/dbscan_3.png) | ![`em_3`](images/em_3.png) | ![`fuzzy_3`](images/fuzzy_3.png) | ![`meanshift_3`](images/meanshift_3.png) | ![`spectral_3`](images/spectral_3.png) | ![`ncuts_3`](images/ncuts_3.png) |
| cluster4 | ![`kmeans_4`](images/kmeans_4.png) | ![`dbscan_4`](images/dbscan_4.png) | ![`em_4`](images/em_4.png) | ![`fuzzy_4`](images/fuzzy_4.png) | ![`meanshift_4`](images/meanshift_4.png) | ![`spectral_4`](images/spectral_4.png) | ![`ncuts_4`](images/ncuts_4.png) |
| cluster5 | ![`kmeans_5`](images/kmeans_5.png) | ![`dbscan_5`](images/dbscan_5.png) | ![`em_5`](images/em_5.png) | ![`fuzzy_5`](images/fuzzy_5.png) | ![`meanshift_5`](images/meanshift_5.png) | ![`spectral_5`](images/spectral_5.png) | ![`ncuts_5`](images/ncuts_5.png) |
| cluster6 | ![`kmeans_6`](images/kmeans_6.png) | ![`dbscan_6`](images/dbscan_6.png) | ![`em_6`](images/em_6.png) | ![`fuzzy_6`](images/fuzzy_6.png) | ![`meanshift_6`](images/meanshift_6.png) | ![`spectral_6`](images/spectral_6.png) | ![`ncuts_6`](images/ncuts_6.png) |

### Parameters

| Algorithm | Parameters |
|-----------|------------|
| dbscan | `--radius 0.08 --minpts 5` (cluster4: `--radius 0.0202 --minpts 5`)¹ |
| em | `--nclusters <K>` |
| fuzzy | `--nclusters <K>` |
| meanshift | `--bandwidth 0.1 --max_iters 30 --tol 0.01` |
| spectral | `--nclusters <K>` |
| ncuts | `--nclusters <K>` |

¹ cluster4 and cluster5 have blob geometry that defeats density-based
separation: their blobs touch via bridge-point chains whose spacing (~0.0018)
is smaller than the blobs' internal point spacing (~0.0040 for cluster5,
~0.0027 for cluster4). Consequences:

- **dbscan** on cluster5 yields one cluster at any `(radius, minpts)` — that
  image is the algorithm's correct answer for that data. cluster4 chains at
  radius 0.08 too; its table image uses a tuned `--radius 0.0202` which
  recovers the 3 ground-truth blobs.
- **meanshift** on cluster4 merges the two upper blobs at bandwidth 0.1 and
  over-fragments at any smaller bandwidth (4+ modes at bw≤0.08), so no
  parameter set yields 3 clean blobs; its 2-cluster image is genuine output.
  Global-structure methods (spectral, ncuts) and parametric methods (kmeans,
  em, fuzzy) are unaffected.

`<K>` per dataset: cluster0=3, cluster1=7, cluster2=2, cluster3=2, cluster4=3, cluster5=3, cluster6=2.

Reproduce one image with:

    ./src/cluster --i data/cluster0.pts --algo kmeans --nclusters 3 --o out_tmp/kmeans_0.ppm && \
    sips -s format png out_tmp/kmeans_0.ppm --out images/kmeans_0.png

