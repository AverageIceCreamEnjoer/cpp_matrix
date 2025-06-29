#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

using ld = long double;

class Matrix final {
 public:
  Matrix() noexcept
      : rows_(0), cols_(0), matrix_(nullptr) {}  // Default constructor

  explicit Matrix(int rows, int cols) : rows_{rows}, cols_{cols} {
    if (rows_ < 0 || cols_ < 0)
      throw std::length_error("Matrix size must be non-negative");  // exception
    matrix_ = std::make_unique<ld[]>(rows_ * cols_);
  }  // matrix rows*cols size

  Matrix(const Matrix& other) : rows_{other.rows_}, cols_{other.cols_} {
    matrix_ = std::make_unique<ld[]>(rows_ * cols_);
    std::copy(other.matrix_.get(), other.matrix_.get() + rows_ * cols_,
              matrix_.get());
  }  // copy

  Matrix(Matrix&& other) noexcept
      : rows_{other.rows_},
        cols_{other.cols_},
        matrix_(std::move(other.matrix_)) {
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }  // move

  ~Matrix() noexcept { Free(); }  // Destructor

  int getRows() const noexcept { return rows_; }  // getter/accessor
  int getCols() const noexcept { return cols_; }  // getter/accessor

  void setRow(int row);
  void setCol(int col);
  Matrix Inverse() const;
  Matrix Transpose() const noexcept;
  void SubMatrix(const Matrix& other);
  void SumMatrix(const Matrix& other);
  void MulNumber(const ld num) noexcept;
  void DivNumber(const ld num);
  void MulMatrix(const Matrix& other);
  ld Determinant() const;
  ld Scalar(const Matrix& other) const;
  Matrix CalcComplements() const;
  ld norm2() const noexcept;
  void print() const;
  int Find_Main_Element(int j) const;
  std::vector<ld> row(int index) const;
  ld Cond_InfinityNorm() const;

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
  void SVD(Matrix& U, Matrix& Sigma, Matrix& V, ld Reduction) const;

  Matrix operator-(const Matrix& other) const {
    Matrix tmp(*this);
    tmp.SubMatrix(other);
    return tmp;
  }

  Matrix operator+(const Matrix& other) const {
    Matrix tmp(*this);
    tmp.SumMatrix(other);
    return tmp;
  }

  Matrix operator*(const Matrix& other) const {
    Matrix tmp(*this);
    tmp.MulMatrix(other);
    return tmp;
  }

  Matrix operator*(const ld num) const noexcept {
    Matrix tmp(*this);
    tmp.MulNumber(num);
    return tmp;
  }

  friend Matrix operator*(const ld num, const Matrix& matrix) noexcept {
    Matrix tmp = matrix * num;
    return tmp;
  }

  Matrix operator/(const ld num) const {
    Matrix tmp(*this);
    tmp.DivNumber(num);
    return tmp;
  }

  friend Matrix operator/(const ld num, const Matrix& matrix) {
    Matrix tmp = matrix / num;
    return tmp;
  }

  Matrix& operator/=(ld num) {
    DivNumber(num);
    return *this;
  }

  Matrix& operator=(const Matrix& other) {
    setRow(other.rows_);
    setCol(other.cols_);
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

  ld& operator()(int row, int col) & {
    return const_cast<ld&>(GetElem(row, col));
  }
  ld& operator()(int row, int col) && = delete;
  const ld& operator()(int row, int col) const& { return GetElem(row, col); }
  const ld& operator()(int row, int col) const&& = delete;

 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  std::unique_ptr<ld[]>
      matrix_;  // Pointer to the memory where the matrix is allocated

  // Private methods
  void Free() noexcept {
    rows_ = 0;
    cols_ = 0;
    matrix_ = nullptr;
  }

  const ld& GetElem(int row, int col) const {
    if (!(row < rows_ && col < cols_ && row >= 0 && col >= 0)) {
      std::cout << row << ' ' << col << std::endl;
      throw std::out_of_range("Index is out of range this matrix");
    }
    return matrix_.get()[row * cols_ + col];
  }

  ld Minor(int row, int col) const {
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
};