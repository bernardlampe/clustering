#ifndef __UTILS_H__
#define __UTILS_H__

#include<fstream>

#include "Exception.h"
#include "Image.h"
#include "Matrix.h"
#include "Vec.h"
#include "cmap.h"
#include "types.h"

#define ABS(N) ((N < 0) ? (-N) : (N))

// Gaussian similarity function
template <typename S>
double similarity(const Matrix<S>& pts, u32 i, u32 j, double sigma = 1.0) {
    double dist2 = 0.0;
    for (u32 d = 0; d < pts.cols(); ++d) {
        double diff = static_cast<double>(pts.get(i, d)) - static_cast<double>(pts.get(j, d));
        dist2 += diff * diff;
    }
    return std::exp(-dist2 / (2 * sigma * sigma));
}

Params_t parse_args(int argc, char **argv) {
  Params_t params;

  for(s32 i = 0; i < argc; i++) {
    char *p = argv[i];
    if (strlen(p) > 2 && p[0] == '-' && p[1] == '-' && i < argc - 1) {
      params[p+2] = argv[++i];
    }
  }

  return params;
}

void readPts(const std::string &fname, u32 &rows, u32 &cols, Matrix<float> &pts) {
  std::ifstream ifile;
  u32 r, c, count = 0;

  ifile.open(fname.c_str());
  if (!ifile)
    throw Exception("could not read from file");

  // count file size
  ifile >> rows >> cols;
  while(!ifile.eof()) {
    ifile >> r >> c;
    count++;
  }
  ifile.close();

  // init and read pts
  pts.init(count, 2);
  ifile.open(fname.c_str());

  count = 0;
  if (!ifile)
    throw Exception("could not read from file");
  ifile >> rows >> cols;
  while(!ifile.eof()) {
    ifile >> r >> c;
    pts.set(count, 0, (float)r);
    pts.set(count, 1, (float)c);
    count++;
  }
  ifile.close();
}

void writeLabelImage(
    const std::string &fname, const std::string &ext,
    const u32 rows, const u32 cols, const u32 nclusters,
    const Vec<u8> &labels, const Matrix<float> &pts) {

  Image<RGB_t> labelImage(rows, cols);
  if (ext == "pts") {
    for(u32 i = 0; i < labels.len(); i++) {
      u32 r = pts.get(i, 0);
      u32 c = pts.get(i, 1);
      u32 l = labels[i];
      labelImage.set(r, c, jetMap[l * int(255/nclusters)]);
    }
  } else {
    for(u32 i = 0; i < labels.len(); i++) {
      u32 l = labels[i];
      labelImage.set(i, jetMap[l * int(255/nclusters)]);
    }
  }
  labelImage.writeToFile(fname);
}

#endif // __UTILS_H__
