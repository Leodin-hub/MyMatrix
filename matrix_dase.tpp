#include "matrix.h"

using namespace leo;

template <class T>
Matrix<T>::Matrix() noexcept : rows_(0), columns_(0), depth_(0) {
  InitErrorText_();
  end_.reset(new value_type[1]);
}

template <class T>
Matrix<T>::Matrix(const size_type x, const size_type y, const size_type z) :
rows_(x), columns_(y), depth_(z) {
  InitErrorText_();
  if (!x || !y || !z) throw std::out_of_range(error_text[ZERO]);
  end_.reset(new value_type[1]);
  NewMatrix_();
}

template <class T>
Matrix<T>::Matrix(const size_type x, const size_type y) :
Matrix(1, x, y) { ; }

template <class T>
Matrix<T>::Matrix(const size_type x, const size_type y, const size_type z, const value_type ***value) :
Matrix(x, y, z, TripValue_(x, y, z, value)) { ; }

template <class T>
Matrix<T>::Matrix(const size_type x, const size_type y, const value_type **value) :
Matrix(1, x, y, TripValue_(x, y, value)) { ; }

template <class T>
Matrix<T>::Matrix(const size_type x, const value_type *value) :
Matrix(1, 1, x, TripValue_(x, value)) { ; }

template <class T>
Matrix<T>::Matrix(const value_type value) : Matrix(1, 1, 1) {
  data_[0][0][0] = value;
}

template <class T>
Matrix<T>::Matrix(const std::initializer_list<value_type> value) :
Matrix(std::vector<value_type>(value)) { ; }

template <class T>
Matrix<T>::Matrix(const size_type x, const size_type y, const size_type z, const trip_vector &vector) {
  if (!(vector.empty() || vector[0].empty() || vector[0][0].empty()))
    ConstructMatrix_(x, y, z, VectorToValue_(vector));
}

template <class T>
Matrix<T>::Matrix(const trip_vector &vector) { 
  if (!(vector.empty() || vector[0].empty() || vector[0][0].empty()))
    ConstructMatrix_(vector.size(), vector[0].size(), vector[0][0].size(), VectorToValue_(vector));
}

template <class T>
Matrix<T>::Matrix(const size_type x, const size_type y, const two_vector &vector) :
Matrix(1, x, y, VectorTriple_(vector)) { ; }

template <class T>
Matrix<T>::Matrix(const two_vector &vector) :
Matrix(VectorTriple_(vector)) { ; }

template <class T>
Matrix<T>::Matrix(const size_type x, const std::vector<value_type> &vector) :
Matrix(1, 1, x, VectorTriple_(vector)) { ; }

template <class T>
Matrix<T>::Matrix(const std::vector<value_type> &vector) :
Matrix(VectorTriple_(vector)) { ; }

template <class T>
Matrix<T>::Matrix(const Matrix &other) :
Matrix(other.rows_, other.columns_, other.depth_, other.data_) { ; }

template <class T>
Matrix<T>::Matrix(Matrix &&other) {
  if (*this == other) throw std::out_of_range(error_text[AFIGEL]);
  rows_ = other.rows_;
  columns_ = other.columns_;
  depth_ = other.depth_;
  end_.swap(other.end_);
  data_.swap(other.data_);
}

template <class T>
typename Matrix<T>::Matrix& Matrix<T>::operator=(const Matrix &other) {
    if (other.rows_)
      ConstructMatrix_(other.rows_, other.columns_, other.depth_);
      ValueToData_(other.data_);
    return *this;
  }

template <class T>
typename Matrix<T>::Matrix& Matrix<T>::operator=(Matrix &&other) {
  if (*this == other) throw std::out_of_range(error_text[AFIGEL]);
  if (!rows_)
    RemoveMatrix_();
  rows_ = other.rows_;
  columns_ = other.columns_;
  depth_ = other.depth_;
  end_.swap(other.end_);
  data_.swap(other.data_);
  return *this;
}

template <class T>
typename Matrix<T>::Matrix& Matrix<T>::operator=(const trip_vector &vector) {
  if (!(vector.empty() || vector[0].empty() || vector[0][0].empty()))
    CopyVector_(vector);
  return *this;
}

template <class T>
typename Matrix<T>::Matrix& Matrix<T>::operator=(const two_vector &vector) {
  return *this = VectorTriple_(vector);;
}

template <class T>
typename Matrix<T>::Matrix& Matrix<T>::operator=(const std::vector<value_type> &vector) {
  return *this = VectorTriple_(vector);;
}

template <class T>
Matrix<T>::~Matrix() { ; }

template <class T>
void Matrix<T>::PushBack(const value_type value) {
  if (rows_ > 1 || columns_ > 1) throw std::out_of_range(error_text[SIZE_MATRIX]);
  if (!rows_) {
    rows_ = columns_ = depth_ = 1;
    NewMatrix_(value);
  } else {
    PushValue_(value);
  }
}

template <class T>
void Matrix<T>::PushBack(const value_type *value) {
  if (!rows_) throw std::out_of_range(error_text[EMPTY]);
  if (rows_ > 1) throw std::out_of_range(error_text[SIZE_MATRIX]);
  PushBack_(depth_, value);
}

template <class T>
void Matrix<T>::PushBack(const size_type x, const value_type *value) {
  if (rows_ > 1 || (depth_ && depth_ != x)) throw std::out_of_range(error_text[SIZE_MATRIX]);
  if (!rows_) {
    data_.reset(new two_data[1]);
    rows_ = 1;
    depth_ = x;
  }
  PushData_(TripValue_(depth_, value)[0][0]);
}

template <class T>
void Matrix<T>::PushBack(const value_type **value) {
  if (!rows_) throw std::out_of_range(error_text[EMPTY]);
  PushBack(columns_, depth_, value);
}

template <class T>
void Matrix<T>::PushBack(const size_type x, const size_type y, const value_type **value) {
  if ((columns_ && columns_ != x) || (depth_ && depth_ != y)) throw std::out_of_range(error_text[SIZE_MATRIX]);
  if (!rows_) {
    columns_ = x;
    depth_ = y;
  }
  PushTwo_(TripValue_(columns_, depth_, value)[0]);
}

template <class T>
void Matrix<T>::PushBack(const two_vector &vector) {
  if ((columns_ && vector.size() != columns_) || (depth_ && vector[0].size() != depth_)) throw std::out_of_range(error_text[SIZE_MATRIX]);
  if (!rows_) {
    columns_ = vector.size();
    depth_ = vector[0].size();
  }
  PushTwo_(VectorToValue_(VectorTriple_(vector))[0]);
}

template <class T>
void Matrix<T>::PushBack(const std::vector<value_type> &vector) {
  if (rows_ > 1 || (depth_ && vector.size() != depth_)) throw std::out_of_range(error_text[SIZE_MATRIX]);
  if (!rows_) {
    data_.reset(new two_data[1]);
    rows_ = 1;
    depth_ = vector.size();
  }
  PushData_(VectorToValue_(VectorTriple_(vector))[0][0]);
}

template <class T>
void Matrix<T>::PushBack(const std::initializer_list<value_type> value) {
  PushBack(std::vector<value_type>(value));
}

template <class T>
void Matrix<T>::PushBack(const Matrix<value_type> &other) {
  if (other.rows_ > 1) throw std::out_of_range(error_text[SIZE_MATRIX]);
  if (rows_) {
    if (other.columns_ > 1 && (columns_ != other.columns_ || depth_ != other.depth_))
      throw std::out_of_range(error_text[SIZE_MATRIX]);
    else if (other.depth_ > 1 && depth_ != other.depth_)
      throw std::out_of_range(error_text[SIZE_MATRIX]);
  }
  if (other.columns_ > 1) {
    if (!rows_) {
      columns_ = other.columns_;
      depth_ = other.depth_;
    }
    PushTwo_(other.data_[0]);
  } else if (other.depth_ > 1) {
    if (!rows_) {
      data_.reset(new two_data[1]);
      rows_ = 1;
      depth_ = other.depth_;
    }
    PushData_(other.data_[0][0]);
  } else {
    if (!rows_) {
      rows_ = columns_ = depth_ = 1;
      NewMatrix_(other.data_[0][0][0]);
    } else {
      PushValue_(other.data_[0][0][0]);
    }
  }
}

template <class T>
Matrix<T>::Matrix(const size_type x, const size_type y, const size_type z, const trip_data &value) {
  InitErrorText_();
  if (!x || !y || !z) throw std::out_of_range(error_text[ZERO]);
  end_.reset(new value_type[1]);
  ConstructMatrix_(x, y, z, value);
}

template <class T>
const typename Matrix<T>::two_data& Matrix<T>::operator[](size_type x) const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  if (x >= rows_) throw std::out_of_range(error_text[SIZE_MATRIX]);
  return data_.get()[x];
}

template <class T>
const typename Matrix<T>::size_type* Matrix<T>::Shape() const noexcept {
  return new size_type[3]{rows_, columns_, depth_};
}

template <class T>
const typename Matrix<T>::size_type Matrix<T>::Size() const noexcept {
  return rows_ * columns_ * depth_;
}

template <class T>
const typename Matrix<T>::size_type Matrix<T>::GetRows() const noexcept {
  return rows_;
}

template <class T>
const typename Matrix<T>::size_type Matrix<T>::GetCollumns() const noexcept {
  return columns_;
}

template <class T>
const typename Matrix<T>::size_type Matrix<T>::GetDepth() const noexcept {
  return depth_;
}

template <class T>
void Matrix<T>::Zeros() {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  CompletionMatrix_(0);
}

template <class T>
void Matrix<T>::Zeros(const size_type x, const size_type y, const size_type z) {
  ConstructMatrix_(x, y, z, 0);
}

template <class T>
void Matrix<T>::Zeros(const size_type x, const size_type y) {
  Zeros(1, x, y);
}

template <class T>
void Matrix<T>::Zeros(const size_type x) {
  Zeros(1, 1, x);
}

template <class T>
void Matrix<T>::Ones() {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  CompletionMatrix_(1);
}

template <class T>
void Matrix<T>::Ones(const size_type x, const size_type y, const size_type z) {
  ConstructMatrix_(x, y, z, 1);
}

template <class T>
void Matrix<T>::Ones(const size_type x, const size_type y) {
  Ones(1, x, y);
}

template <class T>
void Matrix<T>::Ones(const size_type x) {
  Ones(1, 1, x);
}

template <class T>
void Matrix<T>::Completion(const size_type x, const size_type y, const size_type z, const value_type dex) {
  ConstructMatrix_(x, y, z, dex);
}

template <class T>
void Matrix<T>::Completion(const size_type x, const size_type y, const value_type dex) {
  Completion(1, x, y);
}

template <class T>
void Matrix<T>::Completion(const size_type x, const value_type dex) {
  Completion(1, 1, x);
}

template <class T>
void Matrix<T>::CompletionToDo(const value_type tos, const value_type dos) {
  value_type dot = dos;
  if (!dot) dot = rows_ * columns_ * depth_;
  value_type t = tos;
  for (size_type i = 0; i < rows_ && t < dot; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_ && t < dot; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_ && t < dot; ++k)
        r[k] = t++;
    }
  }
}

template <class T>
void Matrix<T>::CompletionEmpty(const value_type dex) {
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        if (!r[k]) r[k] = dex;
    }
  }
}

template <class T>
void Matrix<T>::CompletionEmptyToDo() {
  value_type dex = 1;
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k) {
        auto &t = r[k];
        if (t = dex)
          dex += 1;
        else if (t = 0)
          t = dex;
        else
          dex = t + 1;
      }
    }
  }
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Noise(const value_type x) const {
  Matrix result(*this);
  result.NoiseSet(x);
  return result;
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::NoiseDouble() const {
  Matrix result(*this);
  result.NoiseDoubleSet();
  return result;
}

template <class T>
void Matrix<T>::NoiseSet(const value_type x) {
  srand(time(0));
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] = std::fmod(rand(), x + 1);
    }
  }
}

template <class T>
void Matrix<T>::NoiseDoubleSet() {
  srand(time(0));
  for (size_type i = 0; i < rows_; ++i) {
    two_data &r_r = data_[i];
    for (size_type j = 0; j < columns_; ++j) {
      data_type &r = r_r[j];
      for (size_type k = 0; k < depth_; ++k)
        r[k] = (value_type)rand() / (value_type)RAND_MAX;
    }
  }
}

template <class T>
void Matrix<T>::Resize(const size_type x, const size_type y, const size_type z) {
  if (rows_) {
    data_type temp = Lapsha_();
    const size_type xyz = rows_ * columns_ * depth_;
    RemoveMatrix_();
    rows_ = x; columns_ = y; depth_ = z;
    NewMatrix_();
    data_ = LapshaToData_(x, y, z, xyz, temp);
  } else {
    rows_ = x; columns_ = y; depth_ = z;
    NewMatrix_();
  }
}

template <class T>
void Matrix<T>::Resize(const size_type x, const size_type y) {
  Resize(1, x, y);
}

template <class T>
void Matrix<T>::Resize(const size_type x) {
  Resize(1, 1, x);
}

template <class T>
typename Matrix<T>::trip_vector Matrix<T>::ToTripVector() const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  trip_vector result;
  for (int i = 0; i < rows_; ++i) {
    result.push_back(two_vector());
    auto &r_t = result[i];
    two_data &t_t = data_[i];
    for (int j = 0; j < columns_; ++j) {
      r_t.push_back(std::vector<value_type>());
      auto &r = r_t[j];
      data_type &t = t_t[j];
      for (int k = 0; k < depth_; ++k)
        r.push_back(t[k]);
    }
  }
  return result;
}

template <class T>
typename Matrix<T>::two_vector Matrix<T>::ToTwoVector(const size_type x) const {
  return ToTripVector()[x];
}

template <class T>
std::vector<T> Matrix<T>::ToVector(const size_type x, const size_type y) const {
  return ToTripVector()[x][y];
}

template <class T>
typename Matrix<T>::data_type Matrix<T>::Lapsha() const {
  if (rows_ == 0) throw std::out_of_range(error_text[EMPTY]);
  return Lapsha_();
}

template <class T>
typename Matrix<T>::Matrix Matrix<T>::Unique() const {
  Matrix result(*this);
  result.UniqueSet();
  return result;
}

template <class T>
void Matrix<T>::UniqueSet() {
  if (rows_ > 1)
    UniqueTrip_();
  else if (columns_ > 1)
    UniqueTwo_();
  else if (depth_ > 1)
    UniqueValue_();
  else if (rows_ == 0)
    throw std::out_of_range(error_text[EMPTY]);
}