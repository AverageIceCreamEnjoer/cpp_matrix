#pragma once

#include <gtest/gtest.h>

#include "../matrix.h"

class MatrixTest : public ::testing::Test {
 protected:
  Matrix matrix;

  void SetUp() override {
    matrix = Matrix(2, 2);
    matrix(0, 0) = 1;
    matrix(0, 1) = 2;
    matrix(1, 0) = 2;
    matrix(1, 1) = 1;
  }
};

class MatrixTestNoSetUp : public ::testing::Test {
 protected:
  void SetUp() override {}
};