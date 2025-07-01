#pragma once

#include "matrix.h"

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
  return MatrixRowIterator(MatrixRowIterator::it_, MatrixRowIterator::pos_ + n);
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
  return MatrixRowIterator(MatrixRowIterator::it_, MatrixRowIterator::pos_ - n);
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

template <typename T>
Matrix<T>::iterator Matrix<T>::begin() const noexcept {
  return iterator(this->matrix_.get());
}

template <typename T>
Matrix<T>::iterator Matrix<T>::end() const noexcept {
  return iterator(this->matrix_.get() + this->rows_ * this->cols_,
                  this->rows_ * this->cols_);
}