#include "test.h"

// Тест для пустого конструктора
TEST(MatrixTestNoSetUp, EmptyConstructor) {
  Matrix matrix;
  EXPECT_EQ(matrix.GetRows(), 0);
  EXPECT_EQ(matrix.GetCols(), 0);
}

// Тест для стандратного конструктора
TEST(MatrixTestNoSetUp, SuccessSetUp1) {
  Matrix matrix(1, 1);
  EXPECT_EQ(matrix.GetRows(), 1);
  EXPECT_EQ(matrix.GetCols(), 1);
  matrix(0, 0) = 1;
  EXPECT_EQ(matrix(0, 0), 1);
}

// Тест для стандратного конструктора
TEST(MatrixTestNoSetUp, InvalidSetUp1) {
  EXPECT_THROW(Matrix matrix(-1, -1), std::invalid_argument);
}

// Тест для копирующего конструктора
TEST_F(MatrixTest, CopyConstructor1) {
  Matrix matrix2(matrix);
  EXPECT_EQ(matrix.GetRows(), 2);
  EXPECT_EQ(matrix.GetCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix2.GetRows(), 2);
  EXPECT_EQ(matrix2.GetCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
}

// Тест оператора присваивания копированием
TEST_F(MatrixTest, CopyAssignment1) {
  Matrix matrix2;
  matrix2 = matrix;
  EXPECT_EQ(matrix.GetRows(), 2);
  EXPECT_EQ(matrix.GetCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix2.GetRows(), 2);
  EXPECT_EQ(matrix2.GetCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
}

// Тест для конструктора перемещения
TEST_F(MatrixTest, MoveConstructor1) {
  Matrix matrix2(std::move(matrix));
  EXPECT_EQ(matrix2.GetRows(), 2);
  EXPECT_EQ(matrix2.GetCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
  EXPECT_EQ(matrix.GetRows(), 0);
  EXPECT_EQ(matrix.GetCols(), 0);
  EXPECT_THROW(matrix(0, 0), std::out_of_range);
}

// Тест для оператора присваивания перемещением
TEST_F(MatrixTest, MoveAssignment1) {
  Matrix matrix2;
  matrix2 = std::move(matrix);
  EXPECT_EQ(matrix2.GetRows(), 2);
  EXPECT_EQ(matrix2.GetCols(), 2);
  EXPECT_EQ(matrix2(0, 0), 1);
  EXPECT_EQ(matrix.GetRows(), 0);
  EXPECT_EQ(matrix.GetCols(), 0);
  EXPECT_THROW(matrix(0, 0), std::out_of_range);
}

// Тесты SetRow и SetCol
TEST_F(MatrixTest, SetRow1) {
  matrix.SetRow(1);
  EXPECT_EQ(matrix.GetRows(), 1);
  EXPECT_EQ(matrix.GetCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 1), 2);
  EXPECT_THROW(matrix(1, 0), std::out_of_range);
  matrix.SetRow(2);
  EXPECT_EQ(matrix.GetRows(), 2);
  EXPECT_EQ(matrix.GetCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 1), 2);
  EXPECT_EQ(matrix(1, 0), 0);
  EXPECT_EQ(matrix(1, 1), 0);
  EXPECT_THROW(matrix.SetRow(-1), std::invalid_argument);
}

TEST_F(MatrixTest, SetCol1) {
  matrix.SetCol(1);
  EXPECT_EQ(matrix.GetRows(), 2);
  EXPECT_EQ(matrix.GetCols(), 1);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(1, 0), 2);
  EXPECT_THROW(matrix(0, 1), std::out_of_range);
  matrix.SetCol(2);
  EXPECT_EQ(matrix.GetRows(), 2);
  EXPECT_EQ(matrix.GetCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 1), 0);
  EXPECT_EQ(matrix(1, 0), 2);
  EXPECT_EQ(matrix(1, 1), 0);
  EXPECT_THROW(matrix.SetCol(-1), std::invalid_argument);
}

// Тест для функции транспонирования
TEST_F(MatrixTest, Transpose1) {
  Matrix matrix2 = matrix.Transpose();
  EXPECT_EQ(matrix, matrix2);
}

TEST_F(MatrixTest, Transpose2) {
  matrix.SetCol(1);
  matrix(1, 0) = 3;
  Matrix matrix2 = matrix.Transpose();
  EXPECT_NE(matrix, matrix2);
  for (int i = 0; i < matrix.GetRows(); i++)
    for (int j = 0; j < matrix.GetCols(); j++)
      EXPECT_EQ(matrix(i, j), matrix2(j, i));
}

// Тест детерминанта
TEST_F(MatrixTest, Determinant1) { EXPECT_EQ(matrix.Determinant(), -3); }

TEST_F(MatrixTest, Determinant2) {
  matrix.SetRow(1);
  EXPECT_THROW(matrix.Determinant(), std::logic_error);
}

// Тест матрицы алгебраических дополнений
TEST_F(MatrixTest, CalcComplements1) {
  matrix = matrix.CalcComplements();
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 1), -2);
  EXPECT_EQ(matrix(1, 0), -2);
  EXPECT_EQ(matrix(1, 1), 1);
}

TEST_F(MatrixTest, CalcComplements2) {
  matrix.SetRow(1);
  EXPECT_THROW(matrix.CalcComplements(), std::logic_error);
}

// Тест обратной матрицы
TEST_F(MatrixTest, Inverse1) {
  Matrix matrix2 = matrix * matrix.Inverse();
  EXPECT_EQ(matrix2(0, 0), 1);
  EXPECT_EQ(matrix2(0, 1), 0);
  EXPECT_EQ(matrix2(1, 0), 0);
  EXPECT_EQ(matrix2(1, 1), 1);
}

TEST_F(MatrixTest, Inverse2) {
  matrix.SetRow(1);
  EXPECT_THROW(matrix.Inverse(), std::logic_error);
}

TEST_F(MatrixTest, Inverse3) {
  matrix(0, 1) = 1;
  matrix(1, 1) = 2;
  EXPECT_THROW(matrix.Inverse(), std::underflow_error);
}