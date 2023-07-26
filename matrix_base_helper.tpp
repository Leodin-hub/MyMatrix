#include "matrix.h"

using namespace leo;

template <class T>
void Matrix<T>::InitErrorText_() {
  error_text.reset(new std::string[12]{
    "The determinant is 0",
    "The sizes of the matrices are not equal",
    "The matrix is not square",
    "This calculation is only available for a two-dimensional matrix",
    "The transmitted values go beyond the size of the matrix",
    "The columns of the first matrix should equal the columns of the second matrix",
    "The matrix is empty or not initialized",
    "This function cannot be used to reinitialize the matrix, use another overload or another function",
    "You cannot move the matrix into itself",
    "Cannot be initialized with 0 size",
    "It should be a three-dimensional tensor of the form 1xMx3",
    "The size of the vector should be equal to 3"
  });
}

template <class T>
void Matrix<T>::PushValue_(const value_type value) {
  data_type temp = data_[0][0];
  depth_ += 1;
  data_type &r = data_[0][0];
  r.reset(new value_type[depth_]);
  for (size_type i = 0 ; i < depth_ - 1; ++i)
    r[i] = temp[i];
  r[depth_ - 1] = value;
}

template <class T>
void Matrix<T>::PushData_(const data_type &value) {
  two_data temp = data_[0];
  columns_ += 1;
  two_data &r = data_[0];
  r.reset(new data_type[columns_]);
  for (size_type i = 0; i < columns_ - 1; ++i)
    r[i] = temp[i];
  r[columns_ - 1] = value;
}

template <class T>
void Matrix<T>::PushTwo_(const two_data &value) {
  trip_data temp = data_;
  rows_ += 1;
  data_.reset(new two_data[rows_]);
  for (size_type i = 0; i < rows_ - 1; ++i)
    data_[i] = temp[i];
  data_[rows_ - 1] = value;
}

template <class T>
void Matrix<T>::ConstructMatrix_(const size_type x, const size_type y, const size_type z, const trip_data &value) {
    if (!x || !y || !z) throw std::out_of_range(error_text[ZERO]);
  if (rows_) RemoveMatrix_();
  rows_ = x; columns_ = y; depth_ = z;
  NewMatrix_();
  ValueToData_(value);
}

template <class T>
void Matrix<T>::ConstructMatrix_(const size_type x, const size_type y, const size_type z, const value_type value) {
  if (!x || !y || !z) throw std::out_of_range(error_text[ZERO]);
  if (rows_) RemoveMatrix_();
  rows_ = x; columns_ = y; depth_ = z;
  NewMatrix_(value);
}

template <class T>
template <class K> void Matrix<T>::ValueToData_(const K &value) {
  ValueToData_(rows_, columns_, depth_, value, data_);
}

template <class T>
template <class K> void Matrix<T>::ValueToData_(const size_type x, const size_type y, const size_type z, const K &value, trip_data &data) {
  for (size_type i = 0; i < x; ++i) {
    two_data &r_t = data[i];
    auto &r_k = value[i];
    for (size_type j = 0; j < y; ++j) {
      data_type &r = r_t[j];
      auto &t = r_k[j];
      for (size_type k = 0; k < z; ++k)
        r[k] = t[k];
    }
  }
}

template <class T>
void Matrix<T>::NewMatrix_(const value_type dex) {
  data_ = NewMatrix_(rows_, columns_, depth_, dex);
}

template <class T>
typename Matrix<T>::trip_data Matrix<T>::NewMatrix_(const size_type x, const size_type y, const size_type z, const value_type dex) {
  trip_data result(new two_data[x]);
  for (size_type i = 0; i < x; ++i) {
    two_data &r = result[i];
    r.reset(new data_type[y]);
    for (size_type j = 0; j < y; ++j) {
      data_type &t = r[j];
      t.reset(new value_type[z]{dex});
    }
  }
  return std::move(result);
}

template <class T>
void Matrix<T>::RemoveMatrix_() {
  data_.reset();
  rows_ = columns_ = depth_ = 0;
}

template <class T>
void Matrix<T>::CompletionMatrix_(const value_type dex) {
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_t = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_t[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] = dex;
    }
  }
}

template <class T>
void Matrix<T>::CopyVector_(const trip_vector &vector) {
  ConstructMatrix_(vector.size(), vector[0].size(), vector[0][0].size(), VectorToValue_(vector));
}

template <class T>
typename Matrix<T>::trip_vector Matrix<T>::VectorTriple_(const two_vector &vector) {
  return trip_vector{vector};
}

template <class T>
typename Matrix<T>::trip_vector Matrix<T>::VectorTriple_(const std::vector<value_type> &vector) {
  return VectorTriple_(two_vector{vector});
}

template <class T>
typename Matrix<T>::trip_data Matrix<T>::VectorToValue_(const trip_vector &vector) {
  trip_data result = NewMatrix_(vector.size(), vector[0].size(), vector[0][0].size());
  ValueToData_(vector.size(), vector[0].size(), vector[0][0].size(), vector, result);
  return result;
}

template <class T>
typename Matrix<T>::trip_data Matrix<T>::TripValue_(const size_type x, const size_type y, const size_type z, const value_type ***value) {
  trip_data result = NewMatrix_(x, y, z);
  ValueToData_(x, y, z, value, result);
  return result;
}

template <class T>
typename Matrix<T>::trip_data Matrix<T>::TripValue_(const size_type x, const size_type y, const value_type **value) {
  trip_data result = NewMatrix_(1, x, y);
  two_data &r_t = result[0];
  for (int i = 0; i < x; ++i) {
    data_type &r = r_t[i];
    auto &t = value[i];
    for (int j = 0; j < y; ++j)
      r[j] = t[j];
  }
  return result;
}

template <class T>
typename Matrix<T>::trip_data Matrix<T>::TripValue_(const size_type x, const value_type *value) {
  trip_data result = NewMatrix_(1, 1, x);
  data_type &r = result[0][0];
  for (int i = 0; i < x; ++i)
    r[i] = value[i];
  return result;
}

template <class T>
typename Matrix<T>::data_type Matrix<T>::Lapsha_() const {
  data_type result(new value_type[rows_ * columns_ * depth_]);
  for (size_type i = 0, m = 0; i < rows_; ++i) {
    two_data &r_t = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_t[j];
      for (size_type k = 0; k < depth_; ++k)
        result[m++] = r[k];
    }
  }
  return result;
}

template <class T>
typename Matrix<T>::trip_data Matrix<T>::LapshaToData_(const size_type x, const size_type y, const size_type z, const size_type xyz, const data_type &value) const {
  trip_data result = NewMatrix_(x, y, z);
  for (size_type i = 0, m = 0; i < x && m < xyz; ++i) {
    two_data &r_t = result[i];
    for (size_type j = 0; j < y && m < xyz; ++j) {
      data_type &r = r_t[j];
      for (size_type k = 0; k < z && m < xyz; ++k)
        r[k] = value[m++];
    }
  }
  return result;
}

template <class T>
const bool Matrix<T>::IsQuadratish_(const Matrix &other) const {
  return other.columns_ == other.depth_ && other.rows_ == 1;
}

template <class T>
const bool Matrix<T>::IsEqualMull_(const Matrix &other1, const Matrix &other2) const {
  return (other1.rows_ == other2.rows_ && other1.rows_ == 1) && other1.depth_ == other2.columns_;
}

template <class T>
const bool Matrix<T>::IsEqualShape_(const Matrix &other1, const Matrix &other2) const {
  return other1.rows_ == other2.rows_ && other1.columns_ == other2.columns_ && other1.depth_ == other2.depth_;
}

template <class T>
const bool Matrix<T>::IsEqualTwoData_(const size_type x, const size_type y, const two_data &value1, const two_data &value2) const {
  bool result = true;
  for (size_type i = 0; i < x && result; ++i)
    result = IsEqualData_(y, value1[i], value2[i]);
  return result;
}

template <class T>
const bool Matrix<T>::IsEqualData_(const size_type x, const data_type &value1, const data_type &value2) const {
  bool result = true;
  for (size_type i = 0; i < x && result; ++i)
    result = value1[i] == value2[i];
  return result;
}

template <class T>
void Matrix<T>::UniqueTrip_() {
  std::unique_ptr<bool[]> r(new bool[rows_]{false});
  for (size_type i = 0; i < rows_ - 1; ++i)
    for (size_type j = i + 1; j < rows_; ++j)
      if (!r[i] && IsEqualTwoData_(columns_, depth_, data_[i], data_[j]))
          r[j] = true;
  size_type x = 0;
  for (size_type i = 0; i < rows_; x += !r[i++]);
  if (x < rows_) {
    trip_data result(new two_data[x]);
    for (size_type i = 0, j = 0; i < rows_; ++i)
      if (!r[i])
        result[j++] = data_[i];
    rows_ = x;
    data_ = result;
  }
}

template <class T>
void Matrix<T>::UniqueTwo_() {
  std::unique_ptr<bool[]> r(new bool[columns_]{false});
  two_data &r_t = data_[0];
  for (size_type i = 0; i < columns_ - 1; ++i)
    for (size_type j = i + 1; j < columns_; ++j)
      if (!r[i] && IsEqualData_(depth_, r_t[i], r_t[j]))
          r[j] = true;
  size_type x = 0;
  for (size_type i = 0; i < columns_; x += !r[i++]);
  if (x < columns_) {
    trip_data result(new two_data[1]);
    two_data &r_r = result[0];
    r_r.reset(new data_type[x]);
    for (size_type i = 0, j = 0; i < columns_; ++i)
      if (!r[i])
        r_r[j++] = r_t[i];
    columns_ = x;
    data_ = result;
  }
}

template <class T>
void Matrix<T>::UniqueValue_() {
  std::unique_ptr<bool[]> r(new bool[depth_]{false});
  data_type &r_t = data_[0][0];
  for (size_type i = 0; i < depth_ - 1; ++i)
    for (size_type j = i + 1; j < depth_; ++j)
      if (!r[i] && r_t[i] == r_t[j])
          r[j] = true;
  size_type x = 0;
  for (size_type i = 0; i < depth_; x += !r[i++]);
  if (x < depth_) {
    trip_data result = NewMatrix_(1, 1, x);
    data_type &r_r = result[0][0];
    for (size_type i = 0, j = 0; i < depth_; ++i)
      if (!r[i])
        r_r[j++] = r_t[i];
    depth_ = x;
    data_ = result;
  }
}