template <typename T> void Image<T>::convolve(const float *k, const u32 ksize) {
  s32 center = ksize >> 1;
  Image<float> temp(_width, _height);

  // horizontal pass (rows of this)
  for (u32 h = 0; h < _height; h++) {
    for (u32 w = 0; w < _width; w++) {
      float d = 0.0;
      for (s32 c = -center; c <= center; c++) {
        s32 wp = w + c;
        if (wp >= 0 && wp < (s32)_width)
          d += _data[h * _width + wp] * k[center + c];
      }
      temp.set(w, h, d); // temp is transposed
    }
  }

  // vertical pass over the first pass output, stored in temp at (w, h)
  for (u32 h = 0; h < _height; h++) {
    for (u32 w = 0; w < _width; w++) {
      float d = 0.0;
      for (s32 c = -center; c <= center; c++) {
        s32 hp = h + c;
      if (hp >= 0 && hp < (s32)_height)
          d += temp.get(w, hp) * k[center + c];
      }
      _data[h * _width + w] = (T)d; // destination back into this
    }
  }
}

template <typename T>
void Image<T>::convolve(const float *k, const u32 kheight, const u32 kwidth) {
  Image<float> temp(_height, _width);

  // special case if the ksize is even
  s32 kheightp = (kheight % 2) ? kheight : kheight - 1;
  s32 kwidthp = (kwidth % 2) ? kwidth : kwidth - 1;

  // convolve with 2-D kernel
  for (s32 h = 0; h < _height; h++) {
    for (s32 w = 0; w < _width; w++) {
      float d = 0.0; // kernel accumulator

      // loop over kernel
      for (s32 i = 0; i < (s32)kheight; i++) {
        for (s32 j = 0; j < (s32)kwidth; j++) {
          s32 hp = h + i - (kheightp >> 1);
          s32 wp = w + j - (kwidthp >> 1);

          if (hp >= 0 && hp < _height && wp >= 0 && wp < _width) {
            d += _data[hp * _width + wp] * k[i * kwidth + j];
          }
        }
      }

      temp._data[h * _width + w] = (T)d; // temp is flipped
    }
  }

  // copy convolved image to this image
  *this = temp;
}

template <> Matrix<float> Image<u8>::createpts() {
  // one point per row and cols are num dimensions +2 for spatial position
  Matrix<float> pts(_height * _width, 1);

  for(u32 r = 0, j = 0; r < _height; r++) {
    for(u32 c = 0; c < _width; c++) {
      u32 i = r * _width + c;
      pts[j++] = _data[i];
    }
  }

  return pts;
}

template <> Matrix<float> Image<RGB_t>::createpts() {
  // one point per row and cols are num dimensions +2 for spatial position
  Matrix<float> pts(_height * _width, 3);

  for(u32 r = 0, j = 0; r < _height; r++) {
    for(u32 c = 0; c < _width; c++) {
      u32 i = r * _width + c;
      pts[j++] = _data[i][0];
      pts[j++] = _data[i][1];
      pts[j++] = _data[i][2];
    }
  }

  return pts;
}

template <typename T> Image<T> Image<T>::operator+(const Image<T> &im) const {
  Image<T> temp(_height, _width);

  for (u32 i = 0; i < _height * _width; i++) {
    temp._data[i] = _data[i] + im._data[i];
  }

  return (temp);
}

template <typename T> Image<T> Image<T>::operator-(const Image<T> &im) const {
  Image<T> temp(_height, _width);

  for (u32 i = 0; i < _height * _width; i++) {
    temp._data[i] = _data[i] - im._data[i];
  }

  return (temp);
}

template <typename T> Image<T> Image<T>::operator*(const Image<T> &im) const {
  Image<T> temp(_height, _width);

  for (u32 i = 0; i < _height * _width; i++) {
    temp._data[i] = _data[i] * im._data[i];
  }

  return (temp);
}

template <> void Image<u8>::readFromFile(const std::string &fname) {
  std::ifstream ifile;
  std::string magic, cols, rows, max;
  u8 p;

  ifile.open(fname.c_str());
  if (!ifile) {
    throw Exception("could not read from file");
  }

  ifile >> magic >> cols >> rows >> max;
  if (magic != "P5")
    throw Exception("image needs to be P5 pgm");

  init(atoi(rows.c_str()), atoi(cols.c_str()));
  // after the numeric header the single whitespace byte (the newline) is
  // still in the stream -- consume it, otherwise it becomes pixel [0]
  ifile.get();
  for (u32 i = 0; i < _width * _height; i++) {
    ifile.read((s8 *)&p, sizeof(u8));

    _data[i] = p;
  }
  ifile.close();
}

template <> void Image<RGB_t>::readFromFile(const std::string &fname) {
  std::ifstream ifile;
  std::string magic, cols, rows, max;
  u8 r, g, b;

  ifile.open(fname.c_str());
  if (!ifile) {
    throw Exception("could not read from file");
  }

  ifile >> magic >> cols >> rows >> max;
  if (magic != "P6")
    throw Exception("image needs to be a P6 ppm");
  init(atoi(rows.c_str()), atoi(cols.c_str()));
  // consume the header's trailing newline before raw raster data
  ifile.get();
  for (u32 i = 0; i < _width * _height; i++) {
    ifile.read((s8 *)&r, sizeof(u8));
    ifile.read((s8 *)&g, sizeof(u8));
    ifile.read((s8 *)&b, sizeof(u8));

    _data[i][0] = r;
    _data[i][1] = g;
    _data[i][2] = b;
  }

  ifile.close();
}

/* u8 pixel data must round-trip verbatim for PGM files; the generic writer
 * below contrast-stretches to [0,255] which is right for float signal data
 * but would mutate real 8-bit image values. */
template <> void Image<u8>::writeToFile(const std::string &fname) const {
  std::ofstream ofile;
  u32 numElems = _height * _width;

  if (numElems == 0)
    throw(Exception("cannot write a null image object to file"));

  ofile.open(fname.c_str());
  if (!ofile) {
    throw(Exception("could not open image file for writing"));
  }

  ofile << "P5\n" << _width << " " << _height << "\n255\n";
  for (u32 i = 0; i < numElems; i++) {
    u8 val = (u8)_data[i];
    ofile.write((s8 *)&val, sizeof(u8));
  }
  ofile.close();
}

template <typename T>
void Image<T>::writeToFile(const std::string &fname) const {
  float scaleVal;
  T minVal, maxVal;
  u32 numElems = _height * _width;
  std::ofstream ofile;

  if (numElems == 0) {
    throw(Exception("cannot write a null image object to file"));
  }

  minVal = maxVal = _data[0];
  for (u32 i = 1; i < numElems; i++) {
    if (minVal > _data[i])
      minVal = _data[i];
    if (maxVal < _data[i])
      maxVal = _data[i];
  }

  scaleVal = 255.0 / (minVal < maxVal ? maxVal - minVal : 1.0);

  ofile.open(fname.c_str());
  if (!ofile) {
    throw(Exception("could not open image file for writing"));
  }

  ofile << "P5\n" << _width << " " << _height << "\n255\n";
  for (u32 i = 0; i < numElems; i++) {
    u8 val = (u8)((_data[i] - minVal) * scaleVal + 0.5);
    ofile.write((s8 *)&val, sizeof(u8));
  }

  ofile.close();
}

template <> void Image<RGB_t>::writeToFile(const std::string &fname) const {
  std::ofstream ofile;
  u32 numElems = _height * _width;

  ofile.open(fname.c_str());
  if (!ofile) {
    throw(Exception("unable to open file for writing"));
  }

  ofile << "P6" << std::endl
        << _width << " " << _height << std::endl
        << "255" << std::endl;
  for (u32 i = 0; i < numElems; i++) {
    ofile.write((s8 *)&(_data[i][0]), sizeof(u8));
    ofile.write((s8 *)&(_data[i][1]), sizeof(u8));
    ofile.write((s8 *)&(_data[i][2]), sizeof(u8));
  }

  ofile.close();
}

template <> void Image<Vec2f_t>::writeToFile(const std::string &fname) const {
  std::ofstream ofile;
  s8 spac = 10;

  // allocate space
  u8 *img = new u8[_height * _width];
  memset(img, 0, _height * _width * sizeof(u8));

  // construct the graphical vector field
  for (s32 h = 0; h < (s32)_height; h += spac) {
    for (s32 w = 0; w < (s32)_width; w += spac) {
      s8 ex = w + _data[h * _width + w][0];
      s8 ey = h + _data[h * _width + w][1];

      if (ex >= 0 && ex < (s32)_width && ey >= 0 && ey < (s32)_height) {
        drawLine(w, h, ex, ey, (u8)255, _width, img);
      }
    }
  }

  // open file and check
  ofile.open(fname.c_str());
  if (!ofile) {
    throw(Exception("unable to open file for writing"));
  }

  // write to file
  ofile << "P5" << std::endl
        << _width << " " << _height << std::endl
        << "255" << std::endl;
  ofile.write((s8 *)img, _width * _height * sizeof(u8));
  ofile.close();

  delete[] img;
}
