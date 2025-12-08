template <typename T> double Vec<T>::sum() const {
  double sum = 0;
  for (u32 i = 0; i < _n; i++)
    sum += _data[i];
  return sum;
}

template <typename T> double Vec<T>::prod() const {
  double prod = 0;
  for (u32 i = 0; i < _n; i++)
    prod *= _data[i];
  return prod;
}

template <typename T> double Vec<T>::norm() const {
  double prod = 0;
  for (u32 i = 0; i < _n; i++)
    prod += _data[i] * _data[i];
  return sqrt(prod);
}

template <typename T> Vec<T> Vec<T>::abs() const {
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = std::fabs(_data[i]);
  return v;
}

template <typename T> Vec<T> Vec<T>::operator+(const T &c) const {
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] + c;
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator+=(const T &c) {
  for (u32 i = 0; i < _n; i++)
    _data[i] += c;
  return *this;
}

template <typename T> Vec<T> Vec<T>::operator-(const T &c) const {
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] - c;
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator-=(const T &c) {
  for (u32 i = 0; i < _n; i++)
    _data[i] -= c;
  return *this;
}

template <typename T> Vec<T> Vec<T>::operator*(const T &c) const {
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] * c;
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator*=(const T &c) {
  for (u32 i = 0; i < _n; i++)
    _data[i] *= c;
  return *this;
}

template <typename T> Vec<T> Vec<T>::operator/(const T &c) const {
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] / c;
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator/=(const T &c) {
  for (u32 i = 0; i < _n; i++)
    _data[i] /= c;
  return *this;
}

template <typename T> Vec<T> Vec<T>::operator+(const Vec<T> &m) const {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] + m._data[i];
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator+=(const Vec<T> &m) {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  for (u32 i = 0; i < _n; i++)
    _data[i] += m._data[i];
  return *this;
}

template <typename T> Vec<T> Vec<T>::operator-(const Vec<T> &m) const {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] - m._data[i];
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator-=(const Vec<T> &m) {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  for (u32 i = 0; i < _n; i++)
    _data[i] -= m._data[i];
  return *this;
}

template <typename T> Vec<T> Vec<T>::operator*(const Vec<T> &m) const {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] * m._data[i];
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator*=(const Vec<T> &m) {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  for (u32 i = 0; i < _n; i++)
    _data[i] *= m._data[i];
  return *this;
}

template <typename T> Vec<T> Vec<T>::operator/(const Vec<T> &m) const {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  Vec<T> v(_n);
  for (u32 i = 0; i < _n; i++)
    v._data[i] = _data[i] / m._data[i];
  return v;
}

template <typename T> Vec<T> &Vec<T>::operator/=(const Vec<T> &m) {
  if (_n != m._n)
    throw("cannot compose vectors of different dimension");
  for (u32 i = 0; i < _n; i++)
    _data[i] /= m._data[i];
  return *this;
}
