#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include<string>

#include"Image.h"
#include"Matrix.h"
#include"Vec.h"
#include"types.h"

/*
    - Clustering interface supporting multiple algorithms and up to 255 clusters.
    - Inputs:
        string algo   = "kmeans", "dbscan", "em", "fuzzycmeans", "meanshift",
                        "spectral", "normalizedcuts"
        Matrix<S> pts = the vectors to cluster together, sparse or dense matrix
                        one point per row
    - Outputs:
        u8 labels    = labels 0 to 255, size is number of rows of pts
        u8 nclusters = number of clusters 0 to 255
        Matrix<T>    = vectors representing clusters, one per row same dimension
   as input vectors
*/
template <typename S, typename T>
void cluster(const Params_t &params, const Matrix<S> &pts, Vec<u8> &labels, Matrix<T> &clusters);

/* algorithms implemented one per file */
#include "kmeans.h"
#include "dbscan.h"
#include "em.h"
#include "meanshift.h"
#include "ncuts.h"
/*#include "fuzzy.h"
#include "spectral.h"*/

// one point per row of Matrix pts, supports up to 255 clusters
template <typename S = float, typename T = float>
void cluster(Params_t &params, const Matrix<S> &pts, Vec<u8> &labels, Matrix<T> &clusters) {
  // outer switch for algorithm type
  std::string algo = params["c"];
  if (algo == "kmeans") {

    // check params for kmeans algo
    if (params.find("nclusters") == params.end())
      throw(Exception("kmeans requires --nclusters <int [2 - 255]>"));

    kmeans(pts, std::stoi(params["nclusters"]), labels, clusters);
  }
  else if (algo == "dbscan") {

    // check params for dbscan algo
    if (params.find("radius") == params.end())
      throw(Exception("dbscan requires --radius <float>"));

    if (params.find("minpts") == params.end())
      throw(Exception("dbscan requires --minpts <int>"));

    dbscan(pts, std::stof(params["radius"]), std::stoi(params["minpts"]), labels, clusters);
  }
  else if (algo == "em") {
    if (params.find("nclusters") == params.end())
      throw(Exception("em requires --nclusters <int>"));

    em(pts, std::stoi(params["nclusters"]), labels, clusters);
  }
  else if (algo == "meanshift") {
    if (params.find("bandwidth") == params.end())
      throw(Exception("meanshift requires --bandwidth <double>"));

    if (params.find("max_iters") == params.end())
      throw(Exception("meanshift requires --max_iters <int>"));

    if (params.find("tol") == params.end())
      throw(Exception("meanshift requires --tol <double>"));

    meanshift(pts, std::stod(params["bandwidth"]), std::stoi(params["max_iters"]), std::stod(params["tol"]), labels, clusters);
  }
  else if (algo == "ncuts") {
    if (params.find("nclusters") == params.end())
      throw(Exception("ncuts requires --nclusters <int>"));

    ncuts(pts, std::stoi(params["nclusters"]), labels, clusters);
  }
  else if (algo == "spectral") {
    //spectral(pts, labels, clusters);
  }
  else if (algo == "fuzzy") {
    //fuzzy(pts, labels, clusters);
  }
  else {
    throw(Exception("unknown clustering algorithm"));
  }
}

#endif // __CLUSTER_H__
