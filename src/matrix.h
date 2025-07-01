#pragma once

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

template <typename T>
class Matrix final {
  static_assert(std::is_arithmetic<T>::value,
                "Matrix can only be initialized with arithmetic types");

 public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using size_type = size_t;
  using pointer = T*;

  Matrix() noexcept
      : rows_(0), cols_(0), matrix_(nullptr) {}  // Default constructor

  explicit Matrix(int rows, int cols);  // matrix rows*cols size

  explicit Matrix(
      std::initializer_list<std::initializer_list<value_type>> list);

  Matrix(const Matrix& other);  // copy

  Matrix(Matrix&& other) noexcept;  // move

  ~Matrix() noexcept { Free(); }  // Destructor

  int GetRows() const noexcept { return rows_; }  // getter/accessor
  int GetCols() const noexcept { return cols_; }  // getter/accessor

  void SetRow(int row);
  void SetCol(int col);

 protected:
  template <typename Implementation>
  class MatrixIterator {
   public:
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

   protected:
    pointer it_;
    difference_type pos_;

   public:
    MatrixIterator() = delete;
    MatrixIterator(pointer data, difference_type pos = 0) noexcept
        : it_(data), pos_(pos) {}
    MatrixIterator(const MatrixIterator& other) noexcept
        : it_(other.it_), pos_(other.pos_) {}
    MatrixIterator& operator=(const MatrixIterator& other) noexcept {
      it_ = other.it_;
      pos_ = other.pos_;
      return *this;
    }

    MatrixIterator(MatrixIterator&& other) noexcept
        : it_(other.it_), pos_(other.pos_) {
      other.it_ = nullptr;
      other.pos_ = 0;
    }
    MatrixIterator& operator=(MatrixIterator&& other) noexcept {
      std::swap(it_, other.it_);
      std::swap(pos_, other.pos_);
      return *this;
    }
    ~MatrixIterator() noexcept { it_ = nullptr; };

    Implementation& operator++() { return impl()->operator++(); }
    Implementation operator++(int) { return impl()->operator++(1); }
    Implementation operator+(difference_type n) const {
      return impl()->operator+(n);
    }
    Implementation& operator--() { return impl()->operator--(); }
    Implementation operator--(int) { return impl()->operator--(1); }
    Implementation operator-(difference_type n) const {
      return impl()->operator-(n);
    }
    difference_type operator-(const Implementation& other) const {
      return impl()->operator-(other);
    }
    reference operator*() const { return impl()->operator*(); }
    pointer operator->() const { return impl()->operator->(); }
    bool operator==(const Implementation& other) const {
      return impl()->operator==(other);
    }
    bool operator!=(const Implementation& other) const {
      return impl()->operator!=(other);
    }
    bool operator>(const Implementation& other) const {
      return impl()->operator>(other);
    }
    bool operator<(const Implementation& other) const {
      return impl()->operator<(other);
    }
    bool operator>=(const Implementation& other) const {
      return impl()->operator>=(other);
    }
    bool operator<=(const Implementation& other) const {
      return impl()->operator<=(other);
    }
    reference operator[](difference_type n) const {
      return impl()->operator[](n);
    }

   private:
    Implementation* impl() { return static_cast<Implementation*>(this); }
  };

  class MatrixRowIterator : public MatrixIterator<MatrixRowIterator> {
   public:
    using iterator_category = std::contiguous_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;
    using parent = MatrixIterator<MatrixRowIterator>;

    friend class MatrixIterator<MatrixRowIterator>;
    MatrixRowIterator() = delete;
    MatrixRowIterator(pointer data, difference_type pos = 0) noexcept
        : parent(data, pos) {}
    MatrixRowIterator(const MatrixRowIterator& other) noexcept
        : parent(other.data_, other.pos_) {}
    MatrixRowIterator& operator=(const MatrixRowIterator& other) noexcept {
      parent::it_ = other.it_;
      parent::pos_ = other.pos_;
      return *this;
    }
    MatrixRowIterator(MatrixRowIterator&& other) noexcept {
      parent::it_ = other.it_;
      parent::pos_ = other.pos_;
      other.it_ = nullptr;
      other.pos_ = 0;
    }
    MatrixRowIterator& operator=(MatrixRowIterator&& other) noexcept {
      std::swap(parent::it_, other.it_);
      std::swap(parent::pos_, other.pos_);
      return *this;
    }
    ~MatrixRowIterator() noexcept { parent::it_ = nullptr; };
    MatrixRowIterator& operator++();
    MatrixRowIterator operator++(int);
    MatrixRowIterator operator+(difference_type n) const;
    MatrixRowIterator& operator--();
    MatrixRowIterator operator--(int);
    MatrixRowIterator operator-(difference_type n) const;
    difference_type operator-(const MatrixRowIterator& other) const;
    reference operator*() const;
    pointer operator->() const;
    bool operator==(const MatrixRowIterator& other) const;
    bool operator!=(const MatrixRowIterator& other) const;
    bool operator>(const MatrixRowIterator& other) const;
    bool operator<(const MatrixRowIterator& other) const;
    bool operator>=(const MatrixRowIterator& other) const;
    bool operator<=(const MatrixRowIterator& other) const;
    reference operator[](difference_type n) const;
  };

  class MatrixDiagonalIterator : public MatrixIterator<MatrixDiagonalIterator> {
   public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;
    using parent = MatrixIterator<MatrixDiagonalIterator>;

    size_type cols_ = 0;
    friend class MatrixIterator<MatrixDiagonalIterator>;
    MatrixDiagonalIterator() = delete;
    MatrixDiagonalIterator(pointer data, difference_type pos = 0) noexcept
        : parent(data, pos) {}
    MatrixDiagonalIterator(const MatrixDiagonalIterator& other) noexcept
        : parent(other.data_, other.pos_) {
      cols_ = other.cols_;
    }
    MatrixDiagonalIterator& operator=(
        const MatrixDiagonalIterator& other) noexcept {
      parent::it_ = other.it_;
      parent::pos_ = other.pos_;
      cols_ = other.cols_;
      return *this;
    }
    MatrixDiagonalIterator(MatrixDiagonalIterator&& other) noexcept
        : parent(other) {
      other.it_ = nullptr;
      other.pos_ = 0;
      other.cols_ = 0;
    }
    MatrixDiagonalIterator& operator=(MatrixDiagonalIterator&& other) noexcept {
      std::swap(parent::it_, other.it_);
      std::swap(parent::pos_, other.pos_);
      std::swap(cols_, other.cols_);
      return *this;
    }
    ~MatrixDiagonalIterator() noexcept { parent::it_ = nullptr; };
    MatrixDiagonalIterator& operator++();
    MatrixDiagonalIterator operator++(int);
    MatrixDiagonalIterator operator+(difference_type n) const;
    MatrixDiagonalIterator& operator--();
    MatrixDiagonalIterator operator--(int);
    MatrixDiagonalIterator operator-(difference_type n) const;
    difference_type operator-(const MatrixDiagonalIterator& other) const;
    reference operator*() const;
    pointer operator->() const;
    bool operator==(const MatrixDiagonalIterator& other) const;
    bool operator!=(const MatrixDiagonalIterator& other) const;
    bool operator>(const MatrixDiagonalIterator& other) const;
    bool operator<(const MatrixDiagonalIterator& other) const;
    bool operator>=(const MatrixDiagonalIterator& other) const;
    bool operator<=(const MatrixDiagonalIterator& other) const;
    reference operator[](difference_type n) const;
  };

 public:
  using iterator = MatrixRowIterator;
  iterator begin() const noexcept;
  iterator end() const noexcept;
  MatrixDiagonalIterator d_begin() const;
  MatrixDiagonalIterator d_end() const;
  Matrix Inverse() const;
  Matrix Transpose() const noexcept;
  void SubMatrix(const Matrix& other);
  void SumMatrix(const Matrix& other);
  void MulNumber(const value_type num) noexcept;
  void DivNumber(const value_type num);
  void MulMatrix(const Matrix& other);
  value_type Determinant() const;
  value_type Scalar(const Matrix& other) const;
  Matrix CalcComplements() const;
  value_type Norm2() const noexcept;
  void Print() const;
  std::vector<value_type> Row(int index) const;
  value_type Cond_InfinityNorm() const;

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
  void SVD(Matrix& U, Matrix& Sigma, Matrix& V, value_type Reduction) const;

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

  Matrix operator*(const value_type num) const noexcept {
    Matrix tmp(*this);
    tmp.MulNumber(num);
    return tmp;
  }

  Matrix& operator*=(const value_type num) noexcept {
    MulNumber(num);
    return *this;
  }

  friend Matrix operator*(const value_type num, const Matrix& matrix) noexcept {
    Matrix tmp = matrix * num;
    return tmp;
  }

  Matrix operator/(const value_type num) const {
    Matrix tmp(*this);
    tmp.DivNumber(num);
    return tmp;
  }

  friend Matrix operator/(const value_type num, const Matrix& matrix) {
    return num * matrix.Inverse();
  }

  Matrix& operator/=(value_type num) {
    DivNumber(num);
    return *this;
  }

  Matrix operator/(const Matrix& other) const {
    return *this * other.Inverse();
  }

  Matrix& operator/=(const Matrix& other) { return *this *= other.Inverse(); }

  bool operator==(const Matrix& other) const {
    bool result = (rows_ != other.rows_ || cols_ != other.cols_) ? false : true;
    if (result)
      for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols_; j++)
          if ((*this)(i, j) != other(i, j)) {
            result = false;
            break;
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

  reference operator()(int row, int col) & {
    return const_cast<reference>(GetElem(row, col));
  }
  reference operator()(int row, int col) && = delete;
  const_reference operator()(int row, int col) const& {
    return GetElem(row, col);
  }
  const_reference operator()(int row, int col) const&& = delete;

 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  std::unique_ptr<value_type[]>
      matrix_;  // Pointer to the memory where the matrix is allocated

  // Private methods
  void Free() noexcept {
    rows_ = 0;
    cols_ = 0;
    matrix_ = nullptr;
  }

  const_reference GetElem(int row, int col) const {
    if (!(row < rows_ && col < cols_ && row >= 0 && col >= 0)) {
      throw std::out_of_range("Index " + std::to_string(row) + ' ' +
                              std::to_string(col) +
                              " is out of range this matrix");
    }
    return matrix_.get()[row * cols_ + col];
  }

  value_type Minor(int row, int col) const {
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
#include "matrix_iters.tpp"