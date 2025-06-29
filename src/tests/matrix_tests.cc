#include "test.h"

// Тест для пустого конструктора
TEST(MatrixTestNoSetUp, EmptyConstructor) {
  Matrix matrix;
  EXPECT_EQ(matrix.getRows(), 0);
  EXPECT_EQ(matrix.getCols(), 0);
}

// Тест для стандратного конструктора
TEST(MatrixTestNoSetUp, SuccessSetUp1) {
  Matrix matrix(1, 1);
  EXPECT_EQ(matrix.getRows(), 1);
  EXPECT_EQ(matrix.getCols(), 1);
  matrix(0, 0) = 1;
  EXPECT_EQ(matrix(0, 0), 1);
}

// Тест для стандратного конструктора
TEST(MatrixTestNoSetUp, InvalidSetUp1) {
  EXPECT_THROW(Matrix matrix(-1, -1), std::length_error);
}

// Тест для копирующего конструктора
TEST_F(MatrixTest, CopyConstructor1) {
  Matrix matrix2(matrix);
  EXPECT_EQ(matrix.getRows(), 2);
  EXPECT_EQ(matrix.getCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix2.getRows(), 2);
  EXPECT_EQ(matrix2.getCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
}

// Тест оператора присваивания копированием
TEST_F(MatrixTest, CopyAssignment1) {
  Matrix matrix2;
  matrix2 = matrix;
  EXPECT_EQ(matrix.getRows(), 2);
  EXPECT_EQ(matrix.getCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix2.getRows(), 2);
  EXPECT_EQ(matrix2.getCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
}

// Тест для конструктора перемещения
TEST_F(MatrixTest, MoveConstructor1) {
  Matrix matrix2(std::move(matrix));
  EXPECT_EQ(matrix2.getRows(), 2);
  EXPECT_EQ(matrix2.getCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
  EXPECT_EQ(matrix.getRows(), 0);
  EXPECT_EQ(matrix.getCols(), 0);
  EXPECT_THROW(matrix(0, 0), std::out_of_range);
}

// Тест для оператора присваивания перемещением
TEST_F(MatrixTest, MoveAssignment1) {
  Matrix matrix2;
  matrix2 = std::move(matrix);
  EXPECT_EQ(matrix2.getRows(), 2);
  EXPECT_EQ(matrix2.getCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
  EXPECT_EQ(matrix.getRows(), 0);
  EXPECT_EQ(matrix.getCols(), 0);
  EXPECT_THROW(matrix(0, 0), std::out_of_range);
}
