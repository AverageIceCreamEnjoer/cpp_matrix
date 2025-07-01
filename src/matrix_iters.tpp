#pragma once

#include "matrix.h"

template <typename T>
Matrix<T>::MatrixIterator& Matrix<T>::MatrixIterator::operator++() {
  ++MatrixIterator::pos_;
  ++MatrixIterator::it_;
  return *this;
}

template <typename T>
Matrix<T>::MatrixIterator Matrix<T>::MatrixIterator::operator++(int) {
  MatrixIterator tmp(*this);
  ++(*this);
  return tmp;
}

template <typename T>
Matrix<T>::MatrixIterator Matrix<T>::MatrixIterator::operator+(
    difference_type n) const {
  return MatrixIterator(MatrixIterator::it_, MatrixIterator::pos_ + n);
}

template <typename T>
Matrix<T>::MatrixIterator& Matrix<T>::MatrixIterator::operator--() {
  --MatrixIterator::pos_;
  --MatrixIterator::it_;
  return *this;
}

template <typename T>
Matrix<T>::MatrixIterator Matrix<T>::MatrixIterator::operator--(int) {
  MatrixIterator tmp(*this);
  --(*this);
  return tmp;
}

template <typename T>
Matrix<T>::MatrixIterator Matrix<T>::MatrixIterator::operator-(
    difference_type n) const {
  return MatrixIterator(MatrixIterator::it_, MatrixIterator::pos_ - n);
}

template <typename T>
typename Matrix<T>::MatrixIterator::difference_type
Matrix<T>::MatrixIterator::operator-(const MatrixIterator& other) const {
  return MatrixIterator::pos_ - other.pos_;
}

template <typename T>
Matrix<T>::reference Matrix<T>::MatrixIterator::operator*() const {
  return *MatrixIterator::it_;
}

template <typename T>
Matrix<T>::pointer Matrix<T>::MatrixIterator::operator->() const {
  return MatrixIterator::it_;
}

template <typename T>
bool Matrix<T>::MatrixIterator::operator==(const MatrixIterator& other) const {
  return MatrixIterator::pos_ == other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixIterator::operator!=(const MatrixIterator& other) const {
  return !(*this == other);
}

template <typename T>
bool Matrix<T>::MatrixIterator::operator>(const MatrixIterator& other) const {
  return MatrixIterator::pos_ > other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixIterator::operator<(const MatrixIterator& other) const {
  return MatrixIterator::pos_ < other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixIterator::operator>=(const MatrixIterator& other) const {
  return MatrixIterator::pos_ >= other.pos_;
}

template <typename T>
bool Matrix<T>::MatrixIterator::operator<=(const MatrixIterator& other) const {
  return MatrixIterator::pos_ <= other.pos_;
}

template <typename T>
Matrix<T>::reference Matrix<T>::MatrixIterator::operator[](
    difference_type n) const {
  return *(MatrixIterator::it_ + n);
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