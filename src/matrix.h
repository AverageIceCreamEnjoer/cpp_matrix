#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

template <typename T>
class Matrix final {
  static_assert(std::is_arithmetic<T>::value,
                "Matrix can only be initialized with arithmetic types");

 public:
  Matrix() noexcept
      : rows_(0), cols_(0), matrix_(nullptr) {}  // Default constructor

  explicit Matrix(int rows, int cols);  // matrix rows*cols size

  explicit Matrix(std::initializer_list<std::initializer_list<T>> list);

  Matrix(const Matrix& other);  // copy

  Matrix(Matrix&& other) noexcept;  // move

  ~Matrix() noexcept { Free(); }  // Destructor

  int GetRows() const noexcept { return rows_; }  // getter/accessor
  int GetCols() const noexcept { return cols_; }  // getter/accessor

  void SetRow(int row);
  void SetCol(int col);

 protected:
  class MatrixIterator {
   public:
    using iterator_category = std::contiguous_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

   private:
    pointer data_;
    difference_type pos_;

   public:
    MatrixIterator() = delete;
    MatrixIterator(pointer data, difference_type pos = 0) noexcept
        : data_(data), pos_(pos) {}
    MatrixIterator(const MatrixIterator& other) noexcept
        : data_(other.data_), pos_(other.pos_) {}
    MatrixIterator(MatrixIterator&& other) noexcept
        : data_(other.data_), pos_(other.pos_) {
      other.data_ = nullptr;
      other.pos_ = 0;
    }
    MatrixIterator& operator=(const MatrixIterator& other) noexcept {
      data_ = other.data_;
      pos_ = other.pos_;
      return *this;
    }
    MatrixIterator& operator=(MatrixIterator&& other) noexcept {
      std::swap(data_, other.data_);
      std::swap(pos_, other.pos_);
      return *this;
    }
    ~MatrixIterator() { data_ = nullptr; };

    reference operator*() const noexcept { return data_[pos_]; }
    pointer operator->() const noexcept { return &data_[pos_]; }
    MatrixIterator& operator++() noexcept {
      ++pos_;
      return *this;
    }
    MatrixIterator operator++(int) noexcept {
      MatrixIterator tmp = *this;
      ++pos_;
      return tmp;
    }
    MatrixIterator operator--() noexcept {
      --pos_;
      return *this;
    }
    MatrixIterator operator--(int) noexcept {
      MatrixIterator tmp = *this;
      --pos_;
      return tmp;
    }
    bool operator==(const MatrixIterator& other) const noexcept {
      return data_ == other.data_ && pos_ == other.pos_;
    }
    bool operator!=(const MatrixIterator& other) const noexcept {
      return !(*this == other);
    }
    bool operator>(const MatrixIterator& other) const noexcept {
      return pos_ > other.pos_;
    }
    bool operator>=(const MatrixIterator& other) const noexcept {
      return pos_ >= other.pos_;
    }
    bool operator<(const MatrixIterator& other) const noexcept {
      return pos_ < other.pos_;
    }
    bool operator<=(const MatrixIterator& other) const noexcept {
      return pos_ <= other.pos_;
    }
    difference_type operator-(const MatrixIterator& other) const noexcept {
      return pos_ - other.pos_;
    }
    MatrixIterator operator+(const difference_type offset) const noexcept {
      MatrixIterator tmp = *this;
      tmp.pos_ += offset;
      return tmp;
    }
    MatrixIterator operator-(const difference_type offset) const noexcept {
      MatrixIterator tmp = *this;
      tmp.pos_ -= offset;
      return tmp;
    }
    reference operator[](const difference_type offset) const noexcept {
      return data_[pos_ + offset];
    }
  };

  class ConstMatrixIterator : public MatrixIterator {
   public:
    using MatrixIterator::MatrixIterator;
    using value_type = const T;
    using pointer = const T*;
    using reference = const T&;

    ConstMatrixIterator(const MatrixIterator& other) noexcept
        : MatrixIterator(other) {}

    reference operator*() const noexcept {
      return const_cast<reference>(MatrixIterator::operator*());
    }
    pointer operator->() const noexcept {
      return const_cast<pointer>(MatrixIterator::operator->());
    }
  };

 public:
  using iterator = MatrixIterator;
  using const_iterator = ConstMatrixIterator;
  iterator begin() const noexcept;
  iterator end() const noexcept;
  const_iterator cbegin() const noexcept;
  const_iterator cend() const noexcept;
  Matrix Inverse() const;
  Matrix Transpose() const noexcept;
  void SubMatrix(const Matrix& other);
  void SumMatrix(const Matrix& other);
  void MulNumber(const T num) noexcept;
  void DivNumber(const T num);
  void MulMatrix(const Matrix& other);
  T Determinant() const;
  T Scalar(const Matrix& other) const;
  Matrix CalcComplements() const;
  T Norm2() const noexcept;
  void Print() const;
  std::vector<T> Row(int index) const;
  T Cond_InfinityNorm() const;

  /**
   * @brief Решение системы линейных алгебраических уравнений Ax = f методом
   * LU-разложения
   * @param f правая часть системы
   * @return решение системы
   *
   * Алгоритм LU-разложения. L - нижняя треугольная матрица, U - верхняя
   * треугольная матрица. L и U - матрицы, полученные из матрицы A.
   * A = L * U
   *
   * X = (L^{-1} * (U^{-1} * f))
   */
  Matrix LUSolver(const Matrix& f) const;

  /**
   * @brief QR-разложение на базе отражений Хаусхолдера
   * @return std::pair<Matrix, Matrix> - Q и R из QR-разложения
   *
   * Алгоритм отражений Хаусхолдера.
   * Q - ортогональная матрица, R - верхняя треугольная матрица.
   * Q и R - матрицы, полученные из матрицы A.
   * A = Q * R
   */
  std::pair<Matrix, Matrix> Householder() const;

  /**
   * @brief   Выполнение SVD разложения матрицы A
   * @details
   *          A = U * Sigma * V^T,
   *          где U, V - ортогональные матрицы;
   *          Sigma - диагональная матрица, содержащая сингулярные числа;
   *          V^T - транспонированная матрица V.
   * @param   U  - матрица левых сингулярных векторов
   * @param   Sigma  - матрица сингулярных чисел
   * @param   V  - матрица правых сингулярных векторов
   * @param   Reduction  - величина, до которой производится редукция
   *                       сингулярных чисел
   */
  void SVD(Matrix& U, Matrix& Sigma, Matrix& V, T Reduction) const;

  Matrix operator-() const noexcept {
    Matrix tmp(*this);
    tmp.MulNumber(-1);
    return tmp;
  }

  Matrix operator-(const Matrix& other) const {
    Matrix tmp(*this);
    tmp.SubMatrix(other);
    return tmp;
  }

  Matrix& operator-=(const Matrix& other) {
    SubMatrix(other);
    return *this;
  }

  Matrix operator+(const Matrix& other) const {
    Matrix tmp(*this);
    tmp.SumMatrix(other);
    return tmp;
  }

  Matrix& operator+=(const Matrix& other) {
    SumMatrix(other);
    return *this;
  }

  Matrix operator*(const Matrix& other) const {
    Matrix tmp(*this);
    tmp.MulMatrix(other);
    return tmp;
  }

  Matrix& operator*=(const Matrix& other) {
    MulMatrix(other);
    return *this;
  }

  Matrix operator*(const T num) const noexcept {
    Matrix tmp(*this);
    tmp.MulNumber(num);
    return tmp;
  }

  Matrix& operator*=(const T num) noexcept {
    MulNumber(num);
    return *this;
  }

  friend Matrix operator*(const T num, const Matrix& matrix) noexcept {
    Matrix tmp = matrix * num;
    return tmp;
  }

  Matrix operator/(const T num) const {
    Matrix tmp(*this);
    tmp.DivNumber(num);
    return tmp;
  }

  friend Matrix operator/(const T num, const Matrix& matrix) {
    return num * matrix.Inverse();
  }

  Matrix& operator/=(T num) {
    DivNumber(num);
    return *this;
  }

  Matrix operator/(const Matrix& other) const {
    return *this * other.Inverse();
  }

  Matrix& operator/=(const Matrix& other) { return *this *= other.Inverse(); }

  bool operator==(const Matrix& other) const noexcept {
    bool result = (rows_ != other.rows_ || cols_ != other.cols_) ? false : true;
    if (result) {
      auto ci = other.cbegin();
      for (auto& i : *this)
        if (i != *(ci++)) {
          result = false;
          break;
        }
    }
    return result;
  }

  bool operator!=(const Matrix& other) const { return !(*this == other); }

  Matrix& operator=(const Matrix& other) {
    SetRow(other.rows_);
    SetCol(other.cols_);
    for (int i = 0; i < rows_; i++)
      for (int j = 0; j < cols_; j++) (*this)(i, j) = other(i, j);
    return *this;
  }

  Matrix& operator=(Matrix&& other) {
    if (this != &other) {
      Free();
      std::swap(rows_, other.rows_);
      std::swap(cols_, other.cols_);
      std::swap(matrix_, other.matrix_);
    }
    return *this;
  }

  T& operator()(int row, int col) & {
    return const_cast<T&>(GetElem(row, col));
  }
  T& operator()(int row, int col) && = delete;
  const T& operator()(int row, int col) const& { return GetElem(row, col); }
  const T& operator()(int row, int col) const&& = delete;

 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  std::unique_ptr<T[]>
      matrix_;  // Pointer to the memory where the matrix is allocated

  // Private methods
  void Free() noexcept {
    rows_ = 0;
    cols_ = 0;
    matrix_ = nullptr;
  }

  const T& GetElem(int row, int col) const {
    if (!(row < rows_ && col < cols_ && row >= 0 && col >= 0)) {
      throw std::out_of_range("Index " + std::to_string(row) + ' ' +
                              std::to_string(col) +
                              " is out of range this matrix");
    }
    return matrix_.get()[row * cols_ + col];
  }

  T Minor(int row, int col) const {
    Matrix temp(rows_ - 1, rows_ - 1);
    int shift_i = 0;
    for (int i = 0; i < rows_ - 1; i++) {
      if (i == row) shift_i = 1;
      int shift_j = 0;
      for (int j = 0; j < cols_ - 1; j++) {
        if (j == col) shift_j = 1;
        temp(i, j) = (*this)(i + shift_i, j + shift_j);
      }
    }
    return temp.Determinant();
  }

  int FindMainElement(int j) const;
};

#include "matrix.tpp"