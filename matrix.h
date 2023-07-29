//////////////////////////////////////////////////////////////
/*Моя реализация 3-мерного тензора.
-> Показывает более высокую скорость работы чем
  стандартный 3-мерный вектор.
-> Имеется базовый функционал для работы с
  тензором, матрицами, вектором и скаляром.
-> Реализованны базовые алгебраические операции как
  с числами: суммирование, вычитание, умножение,
  деление, взятие остатка и возведение в степень.
  Так и с тензорами, матрицами и векторами:
  суммирование, вычитание, матричное умножение и деление.
-> Перегруженны большинство методов для работы как с
  объектами этого класса, так и для работы с векторами и
  стандартными массивами.
-> Так же присутствую алгебраические методы для
  работы с квадратными матрицами: вычисление нормали,
  определителя, минора, матрицы дополнений и обратной
  матрицы, преобразование в треугольную и единичную
  матрицу.
-> Так же есть методы для работы с матрицой
  координат типа Nx3: аффинные преобразования и смещение
  координат на заданное значение.
-> Общий для матриц метод транспонирования.
-> Реализованна перегрузка вывода std::cout.
-> Большинство методов дают строгую гарантию безопасности.
-> Методы могут возвращать исключение вида
  std::out_of_range и std::logic_error.*/
//////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <iostream>
#include <memory>

namespace leo {
template <class T> class Matrix {
  public:
  class Iterator;
  using value_type = T;
  using data_type = std::unique_ptr<value_type[]>;
  using two_data = std::unique_ptr<data_type[]>;
  using trip_data = std::unique_ptr<two_data[]>;
  using size_type = unsigned long;
  using two_vector = std::vector<std::vector<value_type>>;
  using trip_vector = std::vector<two_vector>;

  // matrix_iterator.tpp

  class MatrixIterator : public Iterator {
    public:
    MatrixIterator(const trip_data &data, const size_type rows, const size_type columns, const size_type depth, const size_type row, const size_type col, size_type dep, const data_type &end);
    MatrixIterator(const Iterator &other);
    MatrixIterator(Iterator &&other);
    value_type& operator*();
    const value_type* operator->();
    MatrixIterator operator+(const size_type n);
    friend MatrixIterator operator+(const Iterator &iter, const size_type n) {
      return iter + n;
    }
    MatrixIterator operator-(const size_type n);
    friend MatrixIterator operator-(const Iterator &iter, const size_type n) {
      return iter - n;
    }
  };
  
  class ConstMatrixIterator : public Iterator {
    public:
    ConstMatrixIterator(const trip_data &data, const size_type rows, const size_type columns, const size_type depth, const size_type row, const size_type col, size_type dep, const data_type &end);
    ConstMatrixIterator(const Iterator &other);
    ConstMatrixIterator(Iterator &&other);
    const value_type& operator*();
    const value_type* operator->();
    ConstMatrixIterator operator+(const size_type n);
    friend ConstMatrixIterator operator+(const Iterator &iter, const size_type n) {
      return iter + n;
    }
    ConstMatrixIterator operator-(const size_type n);
    friend ConstMatrixIterator operator-(const Iterator &iter, const size_type n) {
      return iter - n;
    }
  };

  using iterator = MatrixIterator;
  using const_iterator = ConstMatrixIterator;

  // Базовый функционал класса реализация в файле:
  // matrix_base.tpp

  // Базовые конструкторы для работы с сырыми данными

  Matrix() noexcept;
  Matrix(const Matrix &other);
  Matrix(Matrix &&other);
  Matrix(const size_type x, const size_type y, const size_type z);
  Matrix(const size_type x, const size_type y);
  // Собирает скаляр из переданного значения
  Matrix(const value_type value);
  // Собирает вектор из переданных значений
  Matrix(const std::initializer_list<value_type> value);

  // Собирает 3-мерный тензор по переданным размерам и заполняет значением из массива.
  // !!!Важно передавать корректные размеры равные или меньше размерам массива!!!
  Matrix(const size_type x, const size_type y, const size_type z, const value_type ***value);
  // Собирает матрица по переданным размерам и заполняет значением из массива.
  // !!!Важно передавать корректные размеры равные или меньше размерам массива!!!
  Matrix(const size_type x, const size_type y, const value_type **value);
  // Собирает вектор по переданным размерам и заполняет значением из массива.
  // !!!Важно передавать корректные размеры равные или меньше размерам массива!!!
  Matrix(const size_type x, const value_type *value);

  Matrix(const size_type x, const size_type y, const size_type z, const trip_vector &vector);
  Matrix(const trip_vector &vector);
  Matrix(const size_type x, const size_type y, const two_vector &vector);
  Matrix(const two_vector &vector);
  Matrix(const size_type x, const std::vector<value_type> &vector);
  Matrix(const std::vector<value_type> &vector);

  Matrix& operator=(const Matrix &other);
  Matrix& operator=(Matrix &&other);
  Matrix& operator=(const trip_vector &vector);
  Matrix& operator=(const two_vector &vector);
  Matrix& operator=(const std::vector<value_type> &vector);

  ~Matrix();

  void PushBack(const value_type value);
  // !!!Важно размеры массива должны совападь с размерами матрицы!!!
  void PushBack(const value_type *value);
  void PushBack(const size_type x, const value_type *value);
  // !!!Важно размеры массива должны совападь с размерами матрицы!!!
  void PushBack(const value_type **value);
  void PushBack(const size_type x, const size_type y, const value_type **value);
  void PushBack(const two_vector &vector);
  void PushBack(const std::vector<value_type> &vector);
  void PushBack(const std::initializer_list<value_type> value);

  void PushBack(const Matrix<value_type> &other);

  const two_data& operator[](size_type x) const;
  const size_type* Shape() const noexcept;
  const size_type Size() const noexcept;
  const size_type GetRows() const noexcept;
  const size_type GetCollumns() const noexcept;
  const size_type GetDepth() const noexcept;


  void Zeros();
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Zeros(const size_type x, const size_type y, const size_type z);
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Zeros(const size_type x, const size_type y);
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Zeros(const size_type x);

  void Ones();
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Ones(const size_type x, const size_type y, const size_type z);
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Ones(const size_type x, const size_type y);
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Ones(const size_type x);

  void Completion(const value_type dex);
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Completion(const size_type x, const size_type y, const size_type z, const value_type dex);
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Completion(const size_type x, const size_type y, const value_type dex);
  // Позволяет пересобрать объект по новым размерам и заполнить значением
  void Completion(const size_type x, const value_type dex);
  void CompletionToDo(const value_type tos = 0, const value_type dos = 0);
  void CompletionEmpty(const value_type dex);
  void CompletionEmptyToDo();

  Matrix Noise(const value_type x) const;
  Matrix NoiseDouble() const;
  void NoiseSet(const value_type x);
  void NoiseDoubleSet();

  void Resize(const size_type x, const size_type y, const size_type z);
  void Resize(const size_type x, const size_type y);
  void Resize(const size_type x);

  void Clear();

  // Преобразует объект в 3-мерный вектор
  trip_vector ToTripVector() const;
  // Преобразует объект в 2-мерный вектор
  two_vector ToTwoVector(const size_type x = 0) const;
  // Преобразует объект в вектор
  std::vector<value_type> ToVector(const size_type x = 0, const size_type y = 0) const;

  // Возвращает массив размера матрицы всех значений объекта
  // Значения будут зачищенны 
  data_type Lapsha() const;

  // Возвращает новый объект уникальных значений текущего
  // Зависит от типа объекта, тензор, матрица или вектор
  Matrix Unique() const;
  // Оставляет только уникальные значения в текущем объекте
  // Зависит от типа объекта, тензор, матрица или вектор
  void UniqueSet();

  friend std::ostream& operator<<(std::ostream &out, const Matrix &curr) {
    if (curr.rows_) {
      if (curr.rows_ > 1) out << "{";
      for (int i = 0; i < curr.rows_; ++i) {
        if (curr.columns_ > 1) out << "{";
        for (int j = 0; j < curr.columns_; ++j) {
          if (curr.depth_ >= 1) out << "{";
          for (int k = 0; k < curr.depth_; ++k) {
            out << curr.data_[i][j][k];
            if (k < curr.depth_ - 1) out << " ";
          }
          if (curr.depth_ >= 1) out << "}";
          if (j < curr.columns_ - 1) out << ", ";
        }
        if (curr.columns_ > 1) out << "}";
        if (i < curr.rows_ - 1) out << "\n";
      }
      if (curr.rows_ > 1) out << "}";
    } else {
      out << "{{{}}}";
    }
    return out;
  }


// Базовые алгебраические методы реализованны в
// matrix_algebraic.tpp

  Matrix operator-() const;

  Matrix operator+(const Matrix &other) const;
  Matrix operator+(const trip_vector &vector) const;
  Matrix operator+(const two_vector &vector) const;
  Matrix operator+(const std::vector<value_type> &vector) const;
  Matrix operator+(const value_type value) const;
  friend Matrix operator+(const trip_vector &vector, const Matrix &other) {
    return other + vector;
  }
  friend Matrix operator+(const two_vector &vector, const Matrix &other) {
    return other + vector;
  }
  friend Matrix operator+(const std::vector<value_type> &vector, const Matrix &other) {
    return other + vector;
  }
  friend Matrix operator+(const value_type value, const Matrix &other) {
    return other + value;
  }
  void operator+=(const Matrix &other);
  void operator+=(const trip_vector &vector);
  void operator+=(const two_vector &vector);
  void operator+=(const std::vector<value_type> &vector);
  void operator+=(const value_type value);

  Matrix operator-(const Matrix &other) const;
  Matrix operator-(const trip_vector &vector) const;
  Matrix operator-(const two_vector &vector) const;
  Matrix operator-(const std::vector<value_type> &vector) const;
  Matrix operator-(const value_type value) const;
  void operator-=(const Matrix &other);
  void operator-=(const trip_vector &vector);
  void operator-=(const two_vector &vector);
  void operator-=(const std::vector<value_type> &vector);
  void operator-=(const value_type value);

  Matrix operator*(const Matrix &other) const;
  Matrix operator*(const two_vector &vector) const;
  Matrix operator*(const std::vector<value_type> &vector) const;
  Matrix operator*(const value_type value) const;
  friend Matrix operator*(const value_type value, const Matrix &other) {
    return other * value;
  }
  void operator*=(const Matrix &other);
  void operator*=(const two_vector &value);
  void operator*=(const std::vector<value_type> &value);
  void operator*=(const value_type value);

  Matrix operator/(const Matrix &other) const;
  Matrix operator/(const two_vector &value) const;
  Matrix operator/(const std::vector<value_type> &value) const;
  Matrix operator/(const value_type value) const;
  void operator/=(const Matrix &other);
  void operator/=(const two_vector &value);
  void operator/=(const std::vector<value_type> &value);
  void operator/=(const value_type value);

  Matrix operator%(const value_type value) const;
  void operator%=(const value_type value);

  Matrix operator^(const value_type value) const;
  void operator^=(const value_type value);

  const bool operator==(const Matrix &other) const;
  const bool operator==(const trip_vector &value) const;
  const bool operator==(const two_vector &value) const;
  const bool operator==(const std::vector<value_type> &value) const;
  const bool operator!=(const Matrix &other) const;
  const bool operator!=(const trip_vector &value) const;
  const bool operator!=(const two_vector &value) const;
  const bool operator!=(const std::vector<value_type> &value) const;

  const value_type Max() const;
  const value_type Min() const;
  const bool Contains(const value_type value) const;

  Matrix Normalization() const;
  void NormalizationSet();

  // Methods for working with a two-dimensional square matrix

  const value_type Normal() const;
  const value_type Determinant() const;
  const value_type MainDiagMull() const;
  Matrix Minor(const size_type x, const size_type y) const;
  Matrix CalcComplementse() const;
  Matrix Inverse() const;
  Matrix Singular() const;
  void SingularSet();
  Matrix TrianMatrix() const;
  void TrianMatrixSet();

  Matrix AffineTransform(const std::vector<value_type> &vector);
  Matrix AffineTransform(const value_type x = 0, const value_type y = 0, const value_type z = 0);
  Matrix AffineTransform(const value_type angel[]);
  void AffineTransformSet(const std::vector<value_type> &vector);
  void AffineTransformSet(const value_type x = 0, const value_type y = 0, const value_type z = 0);
  void AffineTransformSet(const value_type angel[]);
  Matrix Bias(const std::vector<value_type> &vector);
  Matrix Bias(const value_type x = 0, const value_type y = 0, const value_type z = 0);
  Matrix Bias(const value_type bias[]);
  void BiasSet(const std::vector<value_type> &vector);
  void BiasSet(const value_type x = 0, const value_type y = 0, const value_type z = 0);
  void BiasSet(const value_type bias[]);

  Matrix Transpose();
  void TransposeSet();

  // Methods for working with the three-dimensional tensor NxMx3
  // TODO

  // Matrix iteration methods
  // matrix_iterator.tpp

  iterator begin();
  iterator end();

  private:
  size_type rows_, columns_, depth_;
  data_type end_;
  trip_data data_;

  enum ERROR_CODE {
    DETERMINANT,
    SIZE_NOT_EQUAL,
    NOT_SQUARE,
    DIMENSIONAL,
    SIZE_MATRIX,
    COLUMNS_ROWS,
    EMPTY,
    NOT_EMPTY,
    AFIGEL,
    ZERO,
    TENSOR_SIZE,
    VECTOR_SIZE
  };
  std::unique_ptr<std::string[]> error_text;

  // Auxiliary methods for the basic functionality
  // matrix_base_helper.tpp

  void InitErrorText_();

  Matrix(const size_type x, const size_type y, const size_type z, const trip_data &value);

  void PushValue_(const value_type value);
  void PushData_(const data_type &value);
  void PushTwo_(const two_data &value);
  void NewMatrix_(const value_type dex = value_type());
  trip_data NewMatrix_(const size_type x, const size_type y, const size_type z, const value_type dex = value_type());
  void ConstructMatrix_(const size_type x, const size_type y, const size_type z, const trip_data &value);
  void ConstructMatrix_(const size_type x, const size_type y, const size_type z, const value_type value = value_type());
  void RemoveMatrix_();

  template <class K> void ValueToData_(const K &value);
  template <class K> void ValueToData_(const size_type x, const size_type y, const size_type z, const K &value, trip_data &data);
  void CompletionMatrix_(const value_type dex);

  void CopyVector_(const trip_vector &vector);
  trip_vector VectorTriple_(const two_vector &vector);
  trip_vector VectorTriple_(const std::vector<value_type> &vector);
  trip_data VectorToValue_(const trip_vector &vector);

  trip_data TripValue_(const size_type x, const size_type y, const size_type z, const value_type ***value);
  trip_data TripValue_(const size_type x, const size_type y, const value_type **value);
  trip_data TripValue_(const size_type x, const value_type *value);

  data_type Lapsha_() const;
  trip_data LapshaToData_(const size_type x, const size_type y, const size_type z, const size_type xyz, const data_type &value);
  
  const bool IsQuadratish_(const Matrix &other) const;
  const bool IsEqualMull_(const Matrix &other1, const Matrix &other2) const;
  const bool IsEqualShape_(const Matrix &other1, const Matrix &other2) const;
  const bool IsEqualTwoData_(const size_type x, const size_type y, const two_data &value1, const two_data &value2) const;
  const bool IsEqualData_(const size_type x, const data_type &value1, const data_type &value2) const;

  void UniqueTrip_();
  void UniqueTwo_();
  void UniqueValue_();


  // Algebraic methods
  // matrix_algebraic.tpp

  value_type Determ2X2_(const trip_data &data);
  value_type DetermLU_(const size_type x, const trip_data &data);
  bool TrianMatrixEq_(const trip_data &data);
};

// Iterator Methods
// matrix_iterator.tpp

template <class T> class Matrix<T>::Iterator {
  public:
  Iterator(const trip_data &data, const size_type rows, const size_type columns, const size_type depth, const size_type row, const size_type col, size_type dep, const data_type &end);
  Iterator(const Iterator &other);
  Iterator(Iterator &&other);

  const bool operator==(const Iterator &iter);
  const bool operator!=(const Iterator &iter);
  const bool operator>(const Iterator &iter);
  const bool operator>=(const Iterator &iter);
  const bool operator<(const Iterator &iter);
  const bool operator<=(const Iterator &iter);
  value_type* operator++();
  value_type* operator++(int);
  value_type* operator--();
  value_type* operator--(int);
  value_type* operator+=(const size_type n);
  value_type* operator-=(const size_type n);
  const friend size_type operator-(Iterator iter1, Iterator iter2) {
    return iter1.Size_() - iter2.Size_();
  }

  private:
  value_type* NextValue_();
  value_type* PreviousValue_();
  const size_type Size_();

  protected:
  value_type *iterator_;
  size_type rows_, columns_, depth_;
  size_type row_, col_, dep_;
  const trip_data &data_;
  const data_type &end_;
};

}; // namespace leo

#include "matrix_dase.tpp"
#include "matrix_base_helper.tpp"
#include "matrix_algebraic.tpp"
#include "matrix_iterator.tpp"

#endif // MATRIX_H