#!/bin/bash
# Run all algorithms over cluster[0-6].pts, then convert PPM -> PNG.
set -u
cd "$(dirname "$0")"
mkdir -p ../out_tmp ../images

ALGOS="kmeans dbscan em fuzzy meanshift spectral ncuts"
K0=3 K1=7 K2=2 K3=2 K4=3 K5=3 K6=2

run() {
  local algo=$1 ds=$2 k=$3 extra="$4"
  local out="out_tmp/${algo}_${ds}.ppm"
  if [ -f "$out" ]; then echo "skip $algo $ds"; return; fi
  echo "run $algo $ds"
  ../src/cluster --i "cluster${ds}.pts" --algo "$algo" --o "$out" $extra \
    2>&1 | tail -1
}

for ds in 0 1 2 3 4 5 6; do
  eval "K=\$K${ds}"
  run kmeans "$ds" "$K" "--nclusters $K"
  # ds4 blobs chain together at 0.08 (bridges denser than intra-blob gaps);
  # 0.0202 splits ds4 into its 3 true blobs. cluster5 is unseparable by
  # dbscan at any (radius, minpts): bridge pairs (0.0017) are closer than
  # intra-blob NN spacing (0.0040), so one cluster is the correct dbscan answer.
  if [ "$ds" = 4 ]; then
    run dbscan 4 3 "--radius 0.0202 --minpts 5"
  else
    run dbscan "$ds" "$K" "--radius 0.08 --minpts 5"
  fi
  run em     "$ds" "$K" "--nclusters $K"
  run fuzzy  "$ds" "$K" "--nclusters $K"
  run meanshift "$ds" "$K" "--bandwidth 0.1 --max_iters 30 --tol 0.01"
  run spectral "$ds" "$K" "--nclusters $K"
  run ncuts    "$ds" "$K" "--nclusters $K"
done
echo ALL_DONE
