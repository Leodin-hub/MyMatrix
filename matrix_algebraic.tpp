#include "matrix.h"

using namespace leo;


template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator-() const {
  Matrix result(*this);
  result *= (value_type)-1;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator+(const Matrix &other) const {
  Matrix result(*this);
  result += other;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator+(const trip_vector &vector) const {
  Matrix result(*this);
  result += vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator+(const two_vector &vector) const {
  Matrix result(*this);
  result += vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator+(const std::vector<value_type> &vector) const {
  Matrix result(*this);
  result += vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator+(const value_type value) const {
  Matrix result(*this);
  result += value;
  return result;
}

template <class T>
void Matrix<T>::operator+=(const Matrix &other) {
  if (!rows_ || !other.rows_) throw std::out_of_range(error_text[EMPTY]);
  if (!IsEqualShape_(*this, other)) throw std::out_of_range(error_text[SIZE_NOT_EQUAL]);
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    two_data &t_t = other.data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      data_type &t = t_t[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] += t[k];
    }
  }
}

template <class T>
void Matrix<T>::operator+=(const trip_vector &vector) {
  *this += Matrix(vector);
}

template <class T>
void Matrix<T>::operator+=(const two_vector &vector) {
  *this += Matrix(vector);
}

template <class T>
void Matrix<T>::operator+=(const std::vector<value_type> &vector) {
  *this += Matrix(vector);
}

template <class T>
void Matrix<T>::operator+=(const value_type value) {
  if (!rows_) throw std::out_of_range(error_text[EMPTY]);
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] += value;
    }
  }
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator-(const Matrix &other) const {
  Matrix result(*this);
  result -= other;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator-(const trip_vector &vector) const {
  Matrix result(*this);
  result -= vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator-(const two_vector &vector) const {
  Matrix result(*this);
  result -= vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator-(const std::vector<value_type> &vector) const {
  Matrix result(*this);
  result -= vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator-(const value_type value) const {
  Matrix result(*this);
  result -= value;
  return result;
}

template <class T>
void Matrix<T>::operator-=(const Matrix &other) {
  if (!rows_ || !other.rows_) throw std::out_of_range(error_text[EMPTY]);
  if (!IsEqualShape_(*this, other)) throw std::out_of_range(error_text[SIZE_NOT_EQUAL]);
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    two_data &t_t = other.data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      data_type &t = t_t[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] -= t[k];
    }
  }
}

template <class T>
void Matrix<T>::operator-=(const trip_vector &vector) {
  *this -= Matrix(vector);
}

template <class T>
void Matrix<T>::operator-=(const two_vector &vector) {
  *this -= Matrix(vector);
}

template <class T>
void Matrix<T>::operator-=(const std::vector<value_type> &vector) {
  *this -= Matrix(vector);
}

template <class T>
void Matrix<T>::operator-=(const value_type value) {
  *this += -value;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator*(const Matrix &other) const {
  Matrix result(*this);
  result *= other;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator*(const two_vector &vector) const {
  Matrix result(*this);
  result *= vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator*(const std::vector<value_type> &vector) const {
  Matrix result(*this);
  result *= vector;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator*(const value_type value) const {
  Matrix result(*this);
  result *= value;
  return result;
}

template <class T>
void Matrix<T>::operator*=(const Matrix &other) {
  if (!rows_ || !other.rows_) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::logic_error(error_text[DIMENSIONAL]);
  if (!IsEqualMull_(*this, other)) throw std::logic_error(error_text[COLUMNS_ROWS]);
  trip_data result = NewMatrix_(1, columns_, other.depth_);
  two_data &r_t = data_[0];
  two_data &t_t = other.data_[0];
  two_data &r_r = result[0];
  for (size_type i = 0; i < columns_; ++i) {
    data_type &r = r_r[i];
    data_type &t = r_t[i];
    for (size_type j = 0; j < other.depth_; ++j)
      for (size_type k = 0; k < depth_; ++k)
        r[j] += t[k] * t_t[k][j];
  }
  depth_ = other.depth_;
  data_ = std::move(result);
}

template <class T>
void Matrix<T>::operator*=(const two_vector &value) {
  *this *= Matrix(value);
}

template <class T>
void Matrix<T>::operator*=(const std::vector<value_type> &value) {
  *this *= Matrix(value);
}

template <class T>
void Matrix<T>::operator*=(const value_type value) {
  if (!rows_) throw std::out_of_range(error_text[EMPTY]);
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] *= value;
    }
  }
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator/(const Matrix &other) const {
  Matrix result(*this);
  result /= other;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator/(const two_vector &value) const {
  Matrix result(*this);
  result /= value;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator/(const std::vector<value_type> &value) const {
  Matrix result(*this);
  result /= value;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator/(const value_type value) const {
  Matrix result(*this);
  result /= value;
  return result;
}

template <class T>
void Matrix<T>::operator/=(const Matrix &other) {
  this *= other.Inverse();
}

template <class T>
void Matrix<T>::operator/=(const two_vector &value) {
  this /= Matrix(value);
}

template <class T>
void Matrix<T>::operator/=(const std::vector<value_type> &value) {
  this /= Matrix(value);
}

template <class T>
void Matrix<T>::operator/=(const value_type value) {
  *this *= ((value_type)1 / value);
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator%(const value_type value) const {
  Matrix result(*this);
  result %= value;
  return result;
}

template <class T>
void Matrix<T>::operator%=(const value_type value) {
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] = std::fmod(r[k], value);
    }
  }
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::operator^(const value_type value) const {
  Matrix result(*this);
  result ^= value;
  return result;
}

template <class T>
void Matrix<T>::operator^=(const value_type value) {
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] = pow(r[k], value);
    }
  }
}

template <class T>
const bool Matrix<T>::operator==(const Matrix &other) const {
  bool result = rows_ == other.rows_ && columns_ == other.columns_ && depth_ == other.depth_;
  for (size_type i = 0; i < rows_ && result; ++i) {
    two_data &r_r = data_[i];
    two_data  &t_t = other.data_[i];
    for (size_type j = 0; j < columns_ && result; ++j) {
      data_type &r = r_r[j];
      data_type &t = t_t[j];
      for (size_type k = 0; k < depth_ && result; ++k)
        result = r[k] == t[k];
    }
  }
  return result;
}

template <class T>
const bool Matrix<T>::operator==(const trip_vector &value) const {
  return *this == Matrix(value);
}

template <class T>
const bool Matrix<T>::operator==(const two_vector &value) const {
  return *this == Matrix(value);
}

template <class T>
const bool Matrix<T>::operator==(const std::vector<value_type> &value) const {
  return *this == Matrix(value);
}

template <class T>
const bool Matrix<T>::operator!=(const Matrix &other) const {
  return !(*this == other);
}

template <class T>
const bool Matrix<T>::operator!=(const trip_vector &value) const {
  return *this != Matrix(value);
}

template <class T>
const bool Matrix<T>::operator!=(const two_vector &value) const {
  return *this != Matrix(value);
}

template <class T>
const bool Matrix<T>::operator!=(const std::vector<value_type> &value) const {
  return *this != Matrix(value);
}

template <class T>
const typename Matrix<T>::value_type Matrix<T>::Max() const {
  value_type result = value_type();
  if (rows_) {
    result = data_[0][0][0];
    for (size_type i = 0; i < rows_; ++i) {
      two_data &r_r = data_[i];
      for (size_type j = 0; j < columns_; ++j) {
        data_type &r = r_r[j];
        for (size_type k = 0; k < depth_; ++k)
          if (result < r[k]) result = r[k];
      }
    }
  }
  return result;
}

template <class T>
const typename Matrix<T>::value_type Matrix<T>::Min() const {
  value_type result = value_type();
  if (rows_) {
    result = data_[0][0][0];
    for (size_type i = 0; i < rows_; ++i) {
      two_data &r_r = data_[i];
      for (size_type j = 0; j < columns_; ++j) {
        data_type &r = r_r[j];
        for (size_type k = 0; k < depth_; ++k)
          if (result > r[k]) result = r[k];
      }
    }
  }
  return result;
}

template <class T>
const bool Matrix<T>::Contains(const value_type value) const {
  bool result = false;
  if (rows_) {
    for (size_type i = 0; i < rows_ && !result; ++i) {
      two_data &r_r = data_[i];
      for (size_type j = 0; j < columns_ && !result; ++j) {
        data_type &r = r_r[j];
        for (size_type k = 0; k < depth_ && !result; ++k)
          result = value == r[k];
      }
    }
  }
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Normalization() const {
  Matrix result(*this);
  result.NormalizationSet();
  return result;
}

template <class T>
void Matrix<T>::NormalizationSet() {
  if (rows_) {
    *this /= this->Max();
  }
}

template <class T>
const typename Matrix<T>::value_type Matrix<T>::Normal() const {
  value_type result = 0;
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        result += pow(r[k], 2);
    }
  }
  return sqrt(result);
}

template <class T>
const typename Matrix<T>::value_type Matrix<T>::Determinant() const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  if (!IsQuadratish_(*this)) throw std::logic_error(error_text[NOT_SQUARE]);
  value_type result = data_[0][0][0];
  if (columns_ == 2) {
    result = Determ2X2_(data_);
  } else if (MainDiagMull()) {
    result = DetermLU_(columns_, data_);
  } else if (columns_ != 1) {
    data_type &r = data_[0][0];
    for (size_type i = 0; i < columns_; ++i) {
      Matrix temp = Minor(0, i);
      result += r[i] * temp.Determinant() * std::pow(-1, i);
    }
  }
  return result;
}

template <class T>
typename Matrix<T>::value_type Matrix<T>::Determ2X2_(const trip_data &data) {
  return data[0][0][0] * data[0][1][1] - data[0][1][0] * data[0][0][1];
}

template <class T>
typename Matrix<T>::value_type Matrix<T>::DetermLU_(const size_type x, const trip_data &data) {
  Matrix temp(1, x, x, data);
  temp.TrianMatrixSet();
  return temp.MainDiagMull();
}

template <class T>
bool Matrix<T>::TrianMatrixEq_(const trip_data &data) {
  value_type result = 1;
  two_data &r_r = data[0];
  for (size_type i = 0; i < columns_ && result; ++i) {
    data_type &r = r_r[i];
    for (size_type j = 0; j < columns_ && result; ++j)
      if (i - j < 0) result *= r[j];
  }
  return !!result;
}

template <class T>
const typename Matrix<T>::value_type Matrix<T>::MainDiagMull() const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  if (!IsQuadratish_(*this)) throw std::logic_error(error_text[NOT_SQUARE]);
  value_type result = 1;
  two_data &r_r = data_[0];
  for (size_type i = 0; i < columns_; ++i)
    result *= r_r[i][i];
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Minor(const size_type x, const size_type y) const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  // if (!IsQuadratish_(*this)) throw std::logic_error(error_text[NOT_SQUARE]);
  Matrix result(columns_ - 1, depth_ - 1);
  two_data &r_r = data_[0];
  two_data &t_t = result.data_[0];
  for (size_type i = 0, z = 0; i < columns_; ++i) {
    data_type &r = r_r[i];
    for (size_type j = 0; j < depth_; ++j)
      if (i != x && j != y) {
        t_t[z / (columns_ - 1)][z % (columns_ - 1)] = r[j];
        ++z;
      }
  }
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::CalcComplementse() const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  if (!IsQuadratish_(*this)) throw std::logic_error(error_text[NOT_SQUARE]);
  Matrix result(*this);
  if (result.Determinant()) {
    two_data &t_t = result.data_[0];
    for (size_type i = 0; i < columns_; ++i) {
      data_type &t = t_t[i];
      for (size_type j = 0; j < columns_; ++j)
        t[j] = this->Minor(i, j).Determinant() * std::pow(-1, i + j);
    }
  }
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Inverse() const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  if (!IsQuadratish_(*this)) throw std::logic_error(error_text[NOT_SQUARE]);
  Matrix result(this->CalcComplementse());
  value_type determ = this->Determinant();
  result.TransposeSet();
  result *= (value_type)1 / determ;
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Singular() const {
  Matrix result(*this);
  result.SingularSet();
  return result;
}

template <class T>
void Matrix<T>::SingularSet() {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  two_data &r_r = data_[0];
  for (size_type i = 0; i < columns_; ++i)
    r_r[i][i] = 1;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::TrianMatrix() const {
  Matrix result(*this);
  result.TrianMatrixSet();
  return result;
}

template <class T>
void Matrix<T>::TrianMatrixSet() {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  if (!IsQuadratish_(*this)) throw std::logic_error(error_text[NOT_SQUARE]);
  two_data &r_r = data_[0];
  for (int col = 1; col < columns_ && TrianMatrixEq_(data_); col++)
    for (int i = col; i < columns_; ++i) {
      data_type &r = r_r[i];
      data_type &t = r_r[col - 1];
      double f = r[col - 1] / t[col - 1];
      for (int j = col - 1; j < columns_; ++j)
        r[j] -= t[j] * f;
    }
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Transpose() {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  Matrix result(*this);
  result.TransposeSet();
  return result;
}

template <class T>
void Matrix<T>::TransposeSet() {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[DIMENSIONAL]);
  trip_data temp = data_;
  size_type x = columns_;
  rows_ = 1; columns_ = depth_; depth_ = x;
  NewMatrix_();
  two_data &r_r = data_[0];
  two_data &t_t = temp[0];
  for (size_type i = 0; i < depth_; ++i) {
    data_type &t = t_t[i];
    for (size_type j = 0; j < columns_; ++j)
      r_r[j][i] = t[j];
  }
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::AffineTransform(const std::vector<value_type> &vector) {
  Matrix result(*this);
  result.AffineTransformSet(vector);
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::AffineTransform(const value_type x, const value_type y, const value_type z) {
  Matrix result(*this);
  result.AffineTransformSet(x, y, z);
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::AffineTransform(const value_type angel[]) {
  Matrix result(*this);
  result.AffineTransformSet(angel);
  return result;
}

template <class T>
void Matrix<T>::AffineTransformSet(const std::vector<value_type> &vector) {
  if (vector.size() != 3) throw std::out_of_range(error_text[VECTOR_SIZE]);
  this->AffineTransform(vector.data());
}

template <class T>
void Matrix<T>::AffineTransformSet(const value_type x, const value_type y, const value_type z) {
  value_type t[3]{x, y, z};
  this->AffineTransformSet(t);
}

template <class T>
void Matrix<T>::AffineTransformSet(const value_type angel[]) {
  if (!rows_) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1 || depth_ != 3) throw std::logic_error(error_text[TENSOR_SIZE]);
  size_type j = 0, k = 1;
  two_data &result = data_[0];
  for (size_type z = 0; z < 3; ++z) {
    j += z == 2;
    k += z == 1;
    if (angel[z]) {
      value_type angel_r = angel[z] * (3.14 / 180.0);
      for (size_type i = 0; i < columns_; ++i) {
        data_type &r_t = result[i];
        value_type temp_1 = r_t[j];
        value_type temp_2 = r_t[k];
        r_t[j] = std::cos(angel_r) * temp_1 - std::sin(angel_r) * temp_2;
        r_t[k] = std::sin(angel_r) * temp_1 + std::cos(angel_r) * temp_2;
      }
    }
  }
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Bias(const std::vector<value_type> &vector) {
  Matrix result(*this);
  result.BiasSet(vector);
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Bias(const value_type x, const value_type y, const value_type z) {
  Matrix result(*this);
  result.BiasSet(x, y, z);
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Bias(const value_type bias[]) {
  Matrix result(*this);
  result.BiasSet(bias);
  return result;
}

template <class T>
void Matrix<T>::BiasSet(const std::vector<value_type> &vector) {
  if (vector.size() != 3) throw std::out_of_range(error_text[VECTOR_SIZE]);
  this->BiasSet(vector.data());
}

template <class T>
void Matrix<T>::BiasSet(const value_type x, const value_type y, const value_type z) {
  this->BiasSet({x, y, z});
}

template <class T>
void Matrix<T>::BiasSet(const value_type bias[]) {
  if (!rows_) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1 || depth_ != 3) throw std::logic_error(error_text[TENSOR_SIZE]);
  two_data &r_r = data_[0];
  for (size_type i = 0; i < columns_; ++i) {
    data_type &r = r_r[i];
    r[0] += bias[0];
    r[1] += bias[1];
    r[2] += bias[2];
  }
}