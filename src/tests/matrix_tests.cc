#include "test.h"

// Тест для пустого конструктора
TEST(MatrixTestNoSetUp, EmptyConstructor) {
  Matrix<int> matrix;
  EXPECT_EQ(matrix.GetRows(), 0);
  EXPECT_EQ(matrix.GetCols(), 0);
}

// Тест для стандратного конструктора
TEST(MatrixTestNoSetUp, SuccessSetUp1) {
  Matrix<int> matrix(1, 1);
  EXPECT_EQ(matrix.GetRows(), 1);
  EXPECT_EQ(matrix.GetCols(), 1);
  matrix(0, 0) = 1;
  EXPECT_EQ(matrix(0, 0), 1);
}

// Тест для стандратного конструктора
TEST(MatrixTestNoSetUp, InvalidSetUp1) {
  EXPECT_THROW(Matrix<double> matrix(-1, -1), std::invalid_argument);
}

// Тест для списка инициализации
TEST(MatrixTestNoSetUp, ListInit1) {
  Matrix<int> matrix({{1, 2}, {3, 4}});
  EXPECT_EQ(matrix.GetRows(), 2);
  EXPECT_EQ(matrix.GetCols(), 2);
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 1), 2);
  EXPECT_EQ(matrix(1, 0), 3);
  EXPECT_EQ(matrix(1, 1), 4);
}

TEST(MatrixTestNoSetUp, ListInit2) {
  EXPECT_THROW(Matrix<int> matrix({{1, 2}, {3}}), std::invalid_argument);
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
  Matrix<long double> matrix2;
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
  Matrix<long double> matrix2;
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

// Тест разности матриц
TEST_F(MatrixTest, SubMatrix1) {
  matrix.SubMatrix(matrix);
  EXPECT_EQ(matrix(0, 0), 0);
  EXPECT_EQ(matrix(0, 1), 0);
  EXPECT_EQ(matrix(1, 0), 0);
  EXPECT_EQ(matrix(1, 1), 0);
}

TEST_F(MatrixTest, SubMatrix2) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix.SubMatrix(matrix2), std::invalid_argument);
}

TEST_F(MatrixTest, SubMatrix3) {
  matrix = matrix - matrix;
  EXPECT_EQ(matrix(0, 0), 0);
  EXPECT_EQ(matrix(0, 1), 0);
  EXPECT_EQ(matrix(1, 0), 0);
  EXPECT_EQ(matrix(1, 1), 0);
}

TEST_F(MatrixTest, SubMatrix4) {
  matrix = -matrix;
  EXPECT_EQ(matrix(0, 0), -1);
  EXPECT_EQ(matrix(0, 1), -2);
  EXPECT_EQ(matrix(1, 0), -2);
  EXPECT_EQ(matrix(1, 1), -1);
}

TEST_F(MatrixTest, SubMatrix5) {
  matrix -= matrix;
  EXPECT_EQ(matrix(0, 0), 0);
  EXPECT_EQ(matrix(0, 1), 0);
  EXPECT_EQ(matrix(1, 0), 0);
  EXPECT_EQ(matrix(1, 1), 0);
}

TEST_F(MatrixTest, SubMatrix6) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix - matrix2, std::invalid_argument);
}

TEST_F(MatrixTest, SubMatrix7) {
  Matrix matrix2(matrix);
  matrix2.SetCol(1);
  EXPECT_THROW(matrix -= matrix2, std::invalid_argument);
}

// Тест суммы матриц
TEST_F(MatrixTest, SumMatrix1) {
  matrix.SumMatrix(matrix);
  EXPECT_EQ(matrix(0, 0), 2);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 2);
}

TEST_F(MatrixTest, SumMatrix2) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix.SumMatrix(matrix2), std::invalid_argument);
}

TEST_F(MatrixTest, SumMatrix3) {
  matrix = matrix + matrix;
  EXPECT_EQ(matrix(0, 0), 2);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 2);
}

TEST_F(MatrixTest, SumMatrix4) {
  matrix += matrix;
  EXPECT_EQ(matrix(0, 0), 2);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 2);
}

TEST_F(MatrixTest, SumMatrix5) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix + matrix2, std::invalid_argument);
}

TEST_F(MatrixTest, SumMatrix6) {
  Matrix matrix2(matrix);
  matrix2.SetCol(1);
  EXPECT_THROW(matrix += matrix2, std::invalid_argument);
}

// Тест произведения на число
TEST_F(MatrixTest, MulNumber1) {
  matrix.MulNumber(2);
  EXPECT_EQ(matrix(0, 0), 2);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 2);
}

TEST_F(MatrixTest, MulNumber2) {
  matrix.SetRow(1);
  matrix.MulNumber(0);
  EXPECT_EQ(matrix(0, 0), 0);
  EXPECT_EQ(matrix(0, 1), 0);
}

TEST_F(MatrixTest, MulNumber3) {
  matrix = matrix * 2;
  EXPECT_EQ(matrix(0, 0), 2);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 2);
}

TEST_F(MatrixTest, MulNumber4) {
  matrix *= 2;
  EXPECT_EQ(matrix(0, 0), 2);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 2);
}

TEST_F(MatrixTest, MulNumber5) {
  matrix = 2 * matrix;
  EXPECT_EQ(matrix(0, 0), 2);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 2);
}

// Тест деления на число
TEST_F(MatrixTest, DivNumber1) {
  matrix.DivNumber(2);
  EXPECT_EQ(matrix(0, 0), 0.5);
  EXPECT_EQ(matrix(0, 1), 1);
  EXPECT_EQ(matrix(1, 0), 1);
  EXPECT_EQ(matrix(1, 1), 0.5);
}

TEST_F(MatrixTest, DivNumber2) {
  EXPECT_THROW(matrix.DivNumber(0), std::domain_error);
}

TEST_F(MatrixTest, DivNumber3) {
  matrix = matrix / 2;
  EXPECT_EQ(matrix(0, 0), 0.5);
  EXPECT_EQ(matrix(0, 1), 1);
  EXPECT_EQ(matrix(1, 0), 1);
  EXPECT_EQ(matrix(1, 1), 0.5);
}

TEST_F(MatrixTest, DivNumber4) {
  matrix /= 2;
  EXPECT_EQ(matrix(0, 0), 0.5);
  EXPECT_EQ(matrix(0, 1), 1);
  EXPECT_EQ(matrix(1, 0), 1);
  EXPECT_EQ(matrix(1, 1), 0.5);
}

TEST_F(MatrixTest, DivNumber5) {
  EXPECT_THROW(matrix = matrix / 0, std::domain_error);
}

TEST_F(MatrixTest, DivNumber6) { EXPECT_THROW(matrix /= 0, std::domain_error); }

// Тест произведения матриц
TEST_F(MatrixTest, MulMatrix1) {
  matrix.MulMatrix(matrix);
  EXPECT_EQ(matrix(0, 0), 5);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 5);
}

TEST_F(MatrixTest, MulMatrix2) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix.MulMatrix(matrix2), std::invalid_argument);
}

TEST_F(MatrixTest, MulMatrix3) {
  matrix = matrix * matrix;
  EXPECT_EQ(matrix(0, 0), 5);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 5);
}

TEST_F(MatrixTest, MulMatrix4) {
  matrix *= matrix;
  EXPECT_EQ(matrix(0, 0), 5);
  EXPECT_EQ(matrix(0, 1), 4);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 1), 5);
}

TEST_F(MatrixTest, MulMatrix5) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix * matrix2, std::invalid_argument);
}

TEST_F(MatrixTest, MulMatrix6) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix *= matrix2, std::invalid_argument);
}

// Тест деления на матрицу
TEST_F(MatrixTest, DivMatrix1) {
  matrix = matrix / matrix;
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 1), 0);
  EXPECT_EQ(matrix(1, 0), 0);
  EXPECT_EQ(matrix(1, 1), 1);
}

TEST_F(MatrixTest, DivMatrix2) {
  matrix /= matrix;
  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 1), 0);
  EXPECT_EQ(matrix(1, 0), 0);
  EXPECT_EQ(matrix(1, 1), 1);
}

TEST_F(MatrixTest, DivMatrix3) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix / matrix2, std::logic_error);
}

TEST_F(MatrixTest, DivMatrix4) {
  matrix(0, 0) = 2;
  matrix(0, 1) = 1;
  EXPECT_THROW(matrix /= matrix, std::underflow_error);
}

TEST_F(MatrixTest, DivMatrix5) {
  Matrix matrix2 = 1 / matrix;
  matrix2 *= matrix;
  EXPECT_EQ(matrix2(0, 0), 1);
  EXPECT_EQ(matrix2(0, 1), 0);
  EXPECT_EQ(matrix2(1, 0), 0);
  EXPECT_EQ(matrix2(1, 1), 1);
}

// Тест Евклидовой нормы
TEST_F(MatrixTest, Norm2) {
  EXPECT_EQ(matrix.Norm2(), std::sqrt(static_cast<long double>(10)));
}

// Тест скалярного произведения
TEST_F(MatrixTest, Scalar1) {
  Matrix matrix2(matrix);
  EXPECT_FLOAT_EQ(matrix.Scalar(matrix), std::pow(matrix.Norm2(), 2));
}

TEST_F(MatrixTest, Scalar2) {
  Matrix matrix2(matrix);
  matrix2.SetRow(1);
  EXPECT_THROW(matrix.Scalar(matrix2), std::invalid_argument);
}

TEST_F(MatrixTest, Row1) {
  auto row = matrix.Row(0);
  EXPECT_EQ(row[0], 1);
  EXPECT_EQ(row[1], 2);
}

TEST_F(MatrixTest, Row2) { EXPECT_THROW(matrix.Row(-1), std::out_of_range); }

TEST_F(MatrixTest, Iterator1) {
  auto it = matrix.begin();
  EXPECT_EQ(*(it++), 1);
  EXPECT_EQ(*(it++), 2);
  EXPECT_EQ(*(it++), 2);
  EXPECT_EQ(*(it++), 1);
  EXPECT_EQ(it, matrix.end());
}

TEST_F(MatrixTest, Iterator2) {
  auto it = matrix.end();
  EXPECT_EQ(*(--it), 1);
  EXPECT_EQ(*(--it), 2);
  EXPECT_EQ(*(--it), 2);
  EXPECT_EQ(*(--it), 1);
  EXPECT_EQ(it, matrix.begin());
}

TEST_F(MatrixTest, Iterator3) {
  auto it = matrix.begin();
  EXPECT_EQ(*(it + 2), 2);
  EXPECT_EQ(*(it + 3), 1);
}

TEST_F(MatrixTest, Iterator4) {
  auto it = matrix.end();
  EXPECT_EQ(*(it - 2), 2);
  EXPECT_EQ(*(it - 1), 1);
}

TEST_F(MatrixTest, Iterator5) { EXPECT_EQ(matrix.end() - matrix.begin(), 4); }

TEST_F(MatrixTest, Iterator6) { EXPECT_EQ(matrix.begin()[2], 2); }