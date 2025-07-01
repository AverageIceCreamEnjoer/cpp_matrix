#pragma once

#include "matrix.h"

/*
Row Iterator
!----------------------------------------------------------------------------!
*/

template <typename T>
Matrix<T>::MatrixRowIterator& Matrix<T>::MatrixRowIterator::operator++() {
  ++MatrixRowIterator::pos_;
  ++MatrixRowIterator::it_;
  return *this;
}

template <typename T>
Matrix<T>::MatrixRowIterator Matrix<T>::MatrixRowIterator::operator++(int) {
  MatrixRowIterator tmp(*this);
  ++(*this);
  return tmp;
}

template <typename T>
Matrix<T>::MatrixRowIterator Matrix<T>::MatrixRowIterator::operator+(
    difference_type n) const {
  return MatrixRowIterator(MatrixRowIterator::it_ + n,
                           MatrixRowIterator::pos_ + n);
}

template <typename T>
Matrix<T>::MatrixRowIterator& Matrix<T>::MatrixRowIterator::operator--() {
  --MatrixRowIterator::pos_;
  --MatrixRowIterator::it_;
  return *this;
}

template <typename T>
Matrix<T>::MatrixRowIterator Matrix<T>::MatrixRowIterator::operator--(int) {
  MatrixRowIterator tmp(*this);
  --(*this);
  return tmp;
}

template <typename T>
Matrix<T>::MatrixRowIterator Matrix<T>::MatrixRowIterator::operator-(
    difference_type n) const {
  return MatrixRowIterator(MatrixRowIterator::it_ - n,
                           MatrixRowIterator::pos_ - n);
}

template <typename T>
typename Matrix<T>::MatrixRowIterator::difference_type
Matrix<T>::MatrixRowIterator::operator-(const MatrixRowIterator& other) const {
  return MatrixRowIterator::pos_ - other.pos_;
}

template <typename T>
Matrix<T>::reference Matrix<T>::MatrixRowIterator::operator*() const {
  return *MatrixRowIterator::it_;
}

template <typename T>
Matrix<T>::pointer Matrix<T>::MatrixRowIterator::operator->() const {
  return MatrixRowIterator::it_;
}

template <typename T>
bool Matrix<T>::MatrixRowIterator::operator==(
    const MatrixRowIterator& other) const {
  return MatrixRowIterator::pos_ == other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixRowIterator::operator!=(
    const MatrixRowIterator& other) const {
  return !(*this == other);
}

template <typename T>
bool Matrix<T>::MatrixRowIterator::operator>(
    const MatrixRowIterator& other) const {
  return MatrixRowIterator::pos_ > other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixRowIterator::operator<(
    const MatrixRowIterator& other) const {
  return MatrixRowIterator::pos_ < other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixRowIterator::operator>=(
    const MatrixRowIterator& other) const {
  return MatrixRowIterator::pos_ >= other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixRowIterator::operator<=(
    const MatrixRowIterator& other) const {
  return MatrixRowIterator::pos_ <= other.pos_;
}

template <typename T>
Matrix<T>::reference Matrix<T>::MatrixRowIterator::operator[](
    difference_type n) const {
  return *(MatrixRowIterator::it_ + n);
}

/*
Diagonal Iterator
!----------------------------------------------------------------------------!
*/

template <typename T>
Matrix<T>::MatrixDiagonalIterator&
Matrix<T>::MatrixDiagonalIterator::operator++() {
  ++parent::pos_;
  parent::it_ += cols_ + 1;
  return *this;
}

template <typename T>
Matrix<T>::MatrixDiagonalIterator Matrix<T>::MatrixDiagonalIterator::operator++(
    int) {
  MatrixDiagonalIterator tmp(*this);
  ++(*this);
  return tmp;
}

template <typename T>
Matrix<T>::MatrixDiagonalIterator Matrix<T>::MatrixDiagonalIterator::operator+(
    difference_type n) const {
  return MatrixDiagonalIterator(parent::it_ + n * (cols_ + 1),
                                parent::pos_ + n);
}

template <typename T>
Matrix<T>::MatrixDiagonalIterator&
Matrix<T>::MatrixDiagonalIterator::operator--() {
  --parent::pos_;
  parent::it_ -= cols_ + 1;
  return *this;
}

template <typename T>
Matrix<T>::MatrixDiagonalIterator Matrix<T>::MatrixDiagonalIterator::operator--(
    int) {
  MatrixDiagonalIterator tmp(*this);
  --(*this);
  return tmp;
}

template <typename T>
Matrix<T>::MatrixDiagonalIterator Matrix<T>::MatrixDiagonalIterator::operator-(
    difference_type n) const {
  return MatrixDiagonalIterator(parent::it_ - n * (cols_ + 1),
                                parent::pos_ - n);
}

template <typename T>
typename Matrix<T>::MatrixDiagonalIterator::difference_type
Matrix<T>::MatrixDiagonalIterator::operator-(
    const MatrixDiagonalIterator& other) const {
  return parent::pos_ - other.pos_;
}

template <typename T>
Matrix<T>::reference Matrix<T>::MatrixDiagonalIterator::operator*() const {
  return *parent::it_;
}

template <typename T>
Matrix<T>::pointer Matrix<T>::MatrixDiagonalIterator::operator->() const {
  return parent::it_;
}

template <typename T>
bool Matrix<T>::MatrixDiagonalIterator::operator==(
    const MatrixDiagonalIterator& other) const {
  return parent::pos_ == other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixDiagonalIterator::operator!=(
    const MatrixDiagonalIterator& other) const {
  return !(*this == other);
}

template <typename T>
bool Matrix<T>::MatrixDiagonalIterator::operator>(
    const MatrixDiagonalIterator& other) const {
  return parent::pos_ > other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixDiagonalIterator::operator<(
    const MatrixDiagonalIterator& other) const {
  return parent::pos_ < other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixDiagonalIterator::operator>=(
    const MatrixDiagonalIterator& other) const {
  return parent::pos_ >= other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixDiagonalIterator::operator<=(
    const MatrixDiagonalIterator& other) const {
  return parent::pos_ <= other.pos_;
}

template <typename T>
Matrix<T>::reference Matrix<T>::MatrixDiagonalIterator::operator[](
    difference_type n) const {
  return *(parent::it_ + n * (cols_ + 1));
}

template <typename T>
Matrix<T>::iterator Matrix<T>::begin() const noexcept {
  return iterator(this->matrix_.get());
}

template <typename T>
Matrix<T>::iterator Matrix<T>::end() const noexcept {
  return iterator(this->matrix_.get() + this->rows_ * this->cols_,
                  this->rows_ * this->cols_);
}

template <typename T>
Matrix<T>::MatrixDiagonalIterator Matrix<T>::d_begin() const {
  if (this->rows_ != this->cols_)
    throw std::logic_error("MatrixDiagonalIterator: matrix isn't square");
  auto it = MatrixDiagonalIterator(this->matrix_.get());
  it.cols_ = this->cols_;
  return it;
}

template <typename T>
Matrix<T>::MatrixDiagonalIterator Matrix<T>::d_end() const {
  if (this->rows_ != this->cols_)
    throw std::logic_error("MatrixDiagonalIterator: matrix isn't square");
  auto it = MatrixDiagonalIterator(
      this->matrix_.get() + this->rows_ * this->cols_, this->cols_);
  it.cols_ = this->cols_;
  return it;
}