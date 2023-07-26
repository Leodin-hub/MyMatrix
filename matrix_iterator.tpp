#include "matrix.h"

using namespace leo;

template <class T>
typename Matrix<T>::iterator Matrix<T>::begin() {
  return rows_ ? iterator(data_, rows_, columns_, depth_, 0, 0, 0, end_) : end();
}

template <class T>
typename Matrix<T>::iterator Matrix<T>::end() {
  return iterator(data_, rows_, columns_, depth_, rows_, columns_, depth_, end_);
}

template <class T>
Matrix<T>::Iterator::Iterator(const trip_data &data, const size_type rows, const size_type columns, const size_type depth, const size_type row, const size_type col, size_type dep, const data_type &end) :
  data_(data), row_(row), col_(col), dep_(dep), rows_(rows), columns_(columns), depth_(depth), end_(end) {
  iterator_ = row_ == rows_ ? &(end_[0]) : &data_[row_][col_][dep_];
}

template <class T>
Matrix<T>::Iterator::Iterator(const Iterator &other) :
  data_(other.data_), row_(other.row_), col_(other.col_), dep_(other.dep_), rows_(other.rows_), columns_(other.columns_), depth_(other.depth_), end_(other.end_) { ; }

template <class T>
Matrix<T>::Iterator::Iterator(Iterator &&other) :
  data_(other.data_), row_(other.row_), col_(other.col_), dep_(other.dep_), rows_(other.rows_), columns_(other.columns_), depth_(other.depth_), end_(other.end_) { ; }

template <class T>
const bool Matrix<T>::Iterator::operator==(const Iterator &iter) {
  return iterator_ == iter.iterator_;
}

template <class T>
const bool Matrix<T>::Iterator::operator!=(const Iterator &iter) {
  return !(*this == iter);
}

template <class T>
const bool Matrix<T>::Iterator::operator>(const Iterator &iter) {
  if (end_ != iter.end_) std::out_of_range("Sravnivay iterators odnoy matrizi SSUKA!");
  bool result = *this != iter;
  if (!result) {
    if (row_ == rows_) {
      result = false;
    } else if (row_ >= iter.row_) {
      if (row_ == iter.row_) {
        if (col_ >= iter.col_) {
          if (col_ == iter.col_) {
            result = dep_ > iter.dep_;
          } else {
            result = true;
          }
        }
      } else {
        result = true;
      }
    }
  }
  return result;
}

template <class T>
const bool Matrix<T>::Iterator::operator>=(const Iterator &iter) {
  return *this > iter || *this == iter;
}

template <class T>
const bool Matrix<T>::Iterator::operator<(const Iterator &iter) {
  return !(*this > iter);
}

template <class T>
const bool Matrix<T>::Iterator::operator<=(const Iterator &iter) {
  return *this < iter || *this == iter;
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::operator++() {
  return iterator_ = NextValue_();
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::operator++(int) {
  return iterator_ = NextValue_();
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::operator--() {
  return iterator_ = PreviousValue_();
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::operator--(int) {
  return iterator_ = PreviousValue_();
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::operator+=(const size_type n) {
  for (size_type i = 0; i < n; +i) iterator_ = NextValue_();
  return iterator_;
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::operator-=(const size_type n) {
  for (size_type i = 0; i < n; +i) iterator_ = PreviousValue_();
  return iterator_;
}

template <class T>
const typename Matrix<T>::size_type Matrix<T>::Iterator::Size_() {
  size_type result = 0;
  for (size_type i = 0; i <= row_ && i < columns_; ++i)
    for (size_type j = 0; (j <= col_ || (i < row_ && j < columns_)) && j < columns_; ++j)
      for (size_type k = 0; (k <= dep_ || (i < row_ && k < depth_)) && k < depth_; ++k)
        result += 1;
  return result;
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::NextValue_() {
  if (dep_ < depth_ - 1) {
    dep_ += 1;
    iterator_ = &data_[row_][col_][dep_];
  } else if (col_ < columns_ - 1) {
    dep_ = 0;
    col_ += 1;
    iterator_ = &data_[row_][col_][dep_];
  } else if (row_ < rows_ - 1) {
    dep_ = col_ = 0;
    row_ += 1;
    iterator_ = &data_[row_][col_][dep_];
  } else {
    iterator_ = &end_[0];
  }
  return iterator_;
}

template <class T>
typename Matrix<T>::value_type* Matrix<T>::Iterator::PreviousValue_() {
  if (iterator_ == &(end_[0])) {
    row_ = rows_ - 1;
    col_ = columns_ - 1;
    dep_ = depth_ - 1;
  } else if (dep_ > 0) {
    dep_ -= 1;
  } else if (col_ > 0) {
    dep_ = depth_ - 1;
    col_ -= 1;
  } else if (row_ > 0) {
    dep_ = depth_ - 1;
    col_ = columns_ - 1;
    row_ -= 1;
  }
  return iterator_ = &data_[row_][col_][dep_];
}

template <class T>
Matrix<T>::MatrixIterator::MatrixIterator(const trip_data &data, const size_type rows, const size_type columns, const size_type depth, const size_type row, const size_type col, size_type dep, const data_type &end) :
Iterator(data, rows, columns, depth, row, col, dep, end) { ; }

template <class T>
Matrix<T>::MatrixIterator::MatrixIterator(const Iterator &other) :
Iterator(other) { ; }

template <class T>
Matrix<T>::MatrixIterator::MatrixIterator(Iterator &&other) :
Iterator(other) { ; }

template <class T>
typename Matrix<T>::value_type& Matrix<T>::MatrixIterator::operator*() {
  return *this->iterator_;
}

template <class T>
const typename Matrix<T>::value_type* Matrix<T>::MatrixIterator::operator->() {
  return &this->iterator_;
}

template <class T>
typename Matrix<T>::MatrixIterator Matrix<T>::MatrixIterator::operator+(const size_type n) {
  Iterator result(*this);
  result += n;
  return result;
}

template <class T>
typename Matrix<T>::MatrixIterator Matrix<T>::MatrixIterator::operator-(const size_type n) {
  Iterator result(*this);
  result -= n;
  return result;
}

template <class T>
Matrix<T>::ConstMatrixIterator::ConstMatrixIterator(const trip_data &data, const size_type rows, const size_type columns, const size_type depth, const size_type row, const size_type col, size_type dep, const data_type &end) :
Iterator(data, rows, columns, depth, row, col, dep, end) { ; }

template <class T>
Matrix<T>::ConstMatrixIterator::ConstMatrixIterator(const Iterator &other) :
Iterator(other) { ; }

template <class T>
Matrix<T>::ConstMatrixIterator::ConstMatrixIterator(Iterator &&other) :
Iterator(other) { ; }

template <class T>
const typename Matrix<T>::value_type& Matrix<T>::ConstMatrixIterator::operator*() {
  return this->iterator_;
}

template <class T>
const typename Matrix<T>::value_type* Matrix<T>::ConstMatrixIterator::operator->() {
  return &this->iterator_;
}

template <class T>
typename Matrix<T>::ConstMatrixIterator Matrix<T>::ConstMatrixIterator::operator+(const size_type n) {
  Iterator result(*this);
  result += n;
  return result;
}

template <class T>
typename Matrix<T>::ConstMatrixIterator Matrix<T>::ConstMatrixIterator::operator-(const size_type n) {
  Iterator result(*this);
  result -= n;
  return result;
}