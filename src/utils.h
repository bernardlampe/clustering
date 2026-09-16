#ifndef __UTILS_H__
#define __UTILS_H__

#include<fstream>

#include<cstring>
#include<iostream>
#include<cstdlib>
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

// extract-string args; -h/--h is a valueless flag, handle both spellings
// and return early before any value-consumption step can swallow a token
Params_t parse_args(int argc, char **argv) {
  // -h/--h is a valueless flag, handle both spellings and return early
  // before any value-consumption step can swallow the next token
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--h"))
      return Params_t{{"h", ""}};
  }

  Params_t params;

  for(s32 i = 0; i < argc; i++) {
    char *p = argv[i];
    if (strlen(p) > 2 && p[0] == '-' && p[1] == '-' && i < argc - 1) {
      params[p+2] = argv[++i];
    }
  }

  return params;
}


void usage(int exit_code) {
  cout <<
"Usage: cluster --i <in fname [.pgm,.ppm,.pts]> --algo <algo> [algo params] --o <out fname [.ppm]>\n"
"\n"
"    --h            show this help\n"
"    --i <fname>    input file (.pts point list, .ppm and .pgm images)\n"
"    --o <fname>    output colorized label image (.ppm)\n"
"                   (path extension is normalized internally to .ppm)\n"
"\n"
"Algorithms (chosen with --algo; each requires --nclusters unless noted):\n"
"    kmeans, em, fuzzy, spectral, ncuts\n"
"        --nclusters <int [2 - 255]>   number of clusters\n"
"    dbscan\n"
"        --radius <float>   neighborhood radius\n"
"        --minpts <int>     min points to form a dense region\n"
"    meanshift\n"
"        --bandwidth <float>   kernel bandwidth\n"
"        --max_iters <int>     max iterations per point\n"
"        --tol <float>         convergence tolerance\n"
"\n"
"Input format:\n"
"    .pts files start with '<rows> <cols>' then a list of 2D points.\n"
"    .ppm and .pgm images are converted to a point list automatically.\n";
  exit(exit_code);
}

void readPts(const std::string &fname, u32 &rows, u32 &cols, Matrix<float> &pts) {
  std::ifstream ifile;
  u32 r, c, count = 0;

  ifile.open(fname.c_str());
  if (!ifile)
    throw Exception("could not read from file");

  // count file size
  ifile >> rows >> cols;
  // extract-based loop: eof() only trips after a failed read, so the
  // textbook while(!eof){read;count++} pattern overcounts by one
  while (ifile >> r >> c) {
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
  while (ifile >> r >> c) {
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
  // clamp the distinct label count actually used; dbscan labels noise 254
  // and the label->color index must stay inside jetMap[256]
  u32 ncl = nclusters > 255 ? 255 : (nclusters == 0 ? 1 : nclusters);
  u32 stride = 255 / ncl;

  if (ext == "pts") {
    for(u32 i = 0; i < labels.len(); i++) {
      u32 r = pts.get(i, 0);
      u32 c = pts.get(i, 1);
      u32 l = labels[i] < 254 ? labels[i] : 0; // noise painted as cluster 0 edge color
      labelImage.set(r, c, jetMap[l * stride]);
    }
  } else {
    for(u32 i = 0; i < labels.len(); i++) {
      u32 l = labels[i] < 254 ? labels[i] : 0;
      labelImage.set(i, jetMap[l * stride]);
    }
  }
  labelImage.writeToFile(fname);
}

#endif // __UTILS_H__
