#include "matrix.h"

void Matrix::setRow(int row) {
  if (row < 0) throw std::length_error("Matrix size must be non-negative");
  if (row != rows_) {
    Matrix temp(row, cols_);
    int min = std::min(rows_, row);
    for (int i = 0; i < min; i++)
      for (int j = 0; j < cols_; j++) temp(i, j) = (*this)(i, j);
    *this = std::move(temp);
  }
}  // setter/mutator

void Matrix::setCol(int col) {
  if (col < 0) throw std::length_error("Matrix size must be non-negative");
  if (col != cols_) {
    Matrix temp(rows_, col);
    int min = std::min(cols_, col);
    for (int i = 0; i < rows_; i++)
      for (int j = 0; j < min; j++) temp(i, j) = (*this)(i, j);
    *this = std::move(temp);
  }
}  // setter/mutator

Matrix Matrix::Transpose() const noexcept {
  Matrix res(cols_, rows_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) res(j, i) = (*this)(i, j);
  return res;
}

Matrix Matrix::CalcComplements() const {
  if (rows_ != cols_) throw std::out_of_range("Matrix isn't square");
  Matrix res(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res(i, j) = Minor(i, j);
      if ((i + j) % 2) res(i, j) = -res(i, j);
    }
  }
  return res;
}

// Gauss method
ld Matrix::Determinant() const {
  if (rows_ != cols_) throw std::out_of_range("Matrix isn't square");
  ld res = 1;
  Matrix tmp(*this);
  int size = rows_;
  for (int i = 0; i < size; i++) {
    int pivoting = i;
    for (int j = i + 1; j < size; j++)
      if (std::abs(tmp(j, i)) > std::abs(tmp(pivoting, i))) pivoting = j;

    if (std::abs(tmp(pivoting, i)) < 1e-7) {
      res = 0;
      break;
    }
    if (i != pivoting) {
      res = -res;
      for (int j = 0; j < cols_; j++) std::swap(tmp(i, j), tmp(pivoting, j));
    }
    res *= tmp(i, i);
    for (int j = i + 1; j < size; j++) {
      ld koef = tmp(j, i) / tmp(i, i);
      for (int k = i; k < size; k++) tmp(j, k) -= tmp(i, k) * koef;
    }
  }
  return res;
}

Matrix Matrix::Inverse() const {
  if (rows_ != cols_) throw std::out_of_range("Matrix isn't square");
  ld det = Determinant();
  if (std::abs(det) < std::numeric_limits<ld>::min())
    throw std::out_of_range("Determinant must be non-zero");
  return CalcComplements().Transpose() * (1.0 / det);
}

void Matrix::SubMatrix(const Matrix& other) {
  if (!(rows_ == other.rows_ && cols_ == other.cols_))
    throw std::out_of_range("Different matrix sizes");
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) -= other(i, j);
}

void Matrix::SumMatrix(const Matrix& other) {
  if (!(rows_ == other.rows_ && cols_ == other.cols_))
    throw std::out_of_range("Different matrix sizes");
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) += other(i, j);
}

void Matrix::MulNumber(const ld num) noexcept {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) *= num;
}

void Matrix::DivNumber(const ld num) {
  if (std::abs(num) <= std::numeric_limits<ld>::min())
    throw std::out_of_range("Division by zero");
  for (int i = 0; i < rows_; ++i)
    for (int j = 0; j < cols_; ++j) (*this)(i, j) /= num;
}

void Matrix::MulMatrix(const Matrix& other) {
  if (cols_ != other.rows_)
    throw std::out_of_range("Incorrect matrixes size for multiplication");
  Matrix temp(rows_, other.cols_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < other.cols_; j++)
      for (int k = 0; k < cols_; k++) temp(i, j) += (*this)(i, k) * other(k, j);
  *this = std::move(temp);
}

ld Matrix::norm2() const noexcept {
  ld result = 0;
  for (int i = 0; i < rows_; ++i)
    for (int j = 0; j < cols_; ++j) result += std::pow((*this)(i, j), 2);
  result = std::sqrt(result);
  return result;
}

void Matrix::print() const {
  // if (rows_ > 10) throw std::length_error("Matrix is too big.");
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j)
      std::cout << std::fixed << (*this)(i, j) << ' ';
    std::cout << '\n';
  }
  std::cout << std::endl;
}

// поиск ведущего элемента в j-ом столбце матрицы
int Matrix::Find_Main_Element(int j) const {
  int Index = j;
  for (int i = j + 1; i < rows_; i++)
    if (std::abs((*this)(i, j)) > std::abs((*this)(Index, j))) Index = i;
  if (std::abs((*this)(Index, j)) < std::numeric_limits<ld>::min())
    throw std::logic_error("Gauss_Method: degenerate matrix...");
  return Index;
}

ld Matrix::Scalar(const Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::out_of_range("Different matrix sizes");
  ld res = 0;
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) res += (*this)(i, j) * other(i, j);
  return res;
}

std::vector<ld> Matrix::row(int index) const {
  return std::vector<ld>(matrix_.get() + index * cols_,
                         matrix_.get() + (index + 1) * cols_);
}

Matrix Matrix::LUSolver(const Matrix& f) const {
  Matrix LU(*this);
  // инициализация матрицы перестановок строк
  Matrix P(rows_, 1);
  for (int i = 0; i < rows_; i++) P(i, 0) = i;
  // построение верхней треугольной матрицы
  for (int i = 0; i < LU.getRows() - 1; i++) {
    // находим ведущий элемент в i-том столбце
    int I = LU.Find_Main_Element(i);
    // если это не диагональ
    if (I != i) {
      // переставляем строки I и i в СЛАУ местами
      auto Row = LU.row(I);
      for (int j = 0; j < LU.getCols(); ++j) {
        LU(I, j) = LU(i, j);
        LU(i, j) = Row[j];
      }
      int Index = P(i, 0);
      P(i, 0) = P(I, 0);
      P(I, 0) = Index;
    }
    // для оставшихся строк выполним умножение слева на матрицу преобразований
    for (int j = i + 1; j < LU.getRows(); j++) {
      ld help = LU(j, i) / LU(i, i);
      // для уменьшения ошибок вычислений обнуляемые компоненты занулим явно
      LU(j, i) = 0;
      // вычитаем элементы строки i из строк от i + 1 до A.getRows()
      for (int k = i + 1; k < LU.getRows(); k++) {
        LU(j, k) -= help * LU(i, k);
      }
    }
  }
  // построение нижней треугольной матрицы
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < i; j++) {
      ld sum_LikUkj = 0;
      for (int k = 0; k < j; k++) {
        sum_LikUkj += LU(i, k) * LU(k, j);
      }
      LU(i, j) = ((*this)(P(i, 0), j) - sum_LikUkj) / LU(j, j);
    }
  }
  Matrix X(rows_, 1);
  for (int i = 0; i < LU.getRows(); ++i) {
    X(i, 0) = f(P(i, 0), 0);
    for (int j = 0; j < i; ++j) X(i, 0) -= LU(i, j) * X(j, 0);
  }
  for (int i = LU.getRows() - 1; i >= 0; --i) {
    if (std::abs(LU(i, i)) < std::numeric_limits<ld>::min())
      throw std::out_of_range("Division by zero");
    for (int j = i + 1; j < LU.getRows(); ++j) X(i, 0) -= LU(i, j) * X(j, 0);
    X(i, 0) /= LU(i, i);
  }
  return X;
}

std::pair<Matrix, Matrix> Matrix::Householder() const {
  // инициализация вектора отражения
  Matrix p(rows_, 1);
  // вспомогательные переменные
  ld s, beta, mu;
  // запись матрицы А в R
  Matrix R(*this);
  Matrix Q(rows_, cols_);
  for (int i = 0; i < rows_; ++i) Q(i, i) = 1;
  // алгоритм отражений Хаусхолдера
  for (int i = 0; i < R.getCols() - 1; i++) {
    // находим квадрат нормы столбца для обнуления
    s = 0;
    for (int I = i; I < R.getRows(); I++) s += std::pow(R(I, i), 2);
    // если есть ненулевые элементы под диагональю, тогда
    // квадрат нормы вектора для обнуления не совпадает с квадратом
    // диагонального элемента
    if (std::sqrt(std::abs(s - R(i, i) * R(i, i))) >
        std::numeric_limits<ld>::min()) {
      // выбор знака слагаемого beta = sign(-x1)
      if (R(i, i) < 0)
        beta = std::sqrt(s);
      else
        beta = -std::sqrt(s);
      // вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
      mu = 1.0 / (beta * (beta - R(i, i)));
      // формируем вектор p
      for (int I = 0; I < R.getRows(); I++) {
        p(I, 0) = 0;
        if (I >= i) p(I, 0) = R(I, i);
      }
      // изменяем диагональный элемент
      p(i, 0) -= beta;
      // вычисляем новые компоненты матрицы A = Hm * H(m-1) ... * A
      for (int m = i; m < R.getCols(); m++) {
        // произведение S = At * p
        s = 0;
        for (int n = i; n < R.getRows(); n++) {
          s += R(n, m) * p(n, 0);
        }
        s *= mu;
        // A = A - 2 * p * (At * p)^t / ||p||^2
        for (int n = i; n < R.getRows(); n++) {
          R(n, m) -= s * p(n, 0);
        }
      }
      // вычисляем новые компоненты матрицы Q = Q * H1 * H2 * ...
      for (int m = 0; m < R.getRows(); m++) {
        // произведение Q * p
        s = 0;
        for (int n = i; n < R.getRows(); n++) {
          s += Q(m, n) * p(n, 0);
        }
        s *= mu;
        // Q = Q - p * (Q * p)^t
        for (int n = i; n < R.getRows(); n++) {
          Q(m, n) -= s * p(n, 0);
        }
      }
    }
  }
  return std::make_pair(Q, R);
}

ld Matrix::Cond_InfinityNorm() const {
  // проверка на "квадратность" матрицы
  if (rows_ != cols_) throw std::logic_error("Matrix is not square");
  // решатель СЛАУ: A^t = QR и решаем системы A^t * A^(-t) = E
  auto QR = Matrix(*this).Transpose().Householder();
  // проверка на невырожденность
  if (std::abs(QR.second(rows_ - 1, rows_ - 1)) <
      std::numeric_limits<ld>::min())
    throw std::logic_error("Cond(A): detA = 0 ...");
  // максимальные нормы строк (вычисляются на каждом i-ом потоке)
  ld Norma_Row_A = 0;
  ld Norma_Row_A1 = 0;
  // безымянная функция для решения СЛАУ -> столбцы обратной матрицы
  // строка обратной матрицы
  Matrix A1(rows_, 1);
  ld S1, S2;
  // первая и последняя обрабатываемые строки для потока
  int Begin = 0;
  int End = Begin + rows_;
  // решаем системы A^t * A^(-t) = E
  for (int i = Begin; i < End; i++) {
    A1(i, 0) = 1.0;
    Matrix y = QR.first.Transpose() * A1;
    for (int i = A1.getRows() - 1; i >= 0; --i) {
      ld s = 0;
      for (int j = i + 1; j < A1.getRows() - 1; ++j) s += QR.second(i, j);
      A1(i, 0) = (y(i, 0) - s) / QR.second(i, i);
    }
    S1 = 0;
    S2 = 0;
    for (int j = 0; j < rows_; j++) {
      S1 += std::abs((*this)(i, j));
      S2 += std::abs(A1(j, 0));
      A1(j, 0) = 0.0;
    }
    if (Norma_Row_A < S1) Norma_Row_A = S1;
    if (Norma_Row_A1 < S2) Norma_Row_A1 = S2;
  }
  return Norma_Row_A * Norma_Row_A1;
}

void Householder_Col_Transform(Matrix& A, Matrix& U, int i, int j) {
  // вектор отражения
  Matrix p(A.getRows(), 1);
  // вспомогательные переменные
  ld s, beta, mu;
  // находим квадрат нормы столбца для обнуления
  s = 0;
  for (int I = j; I < A.getRows(); I++) s += std::pow(A(I, i), 2);
  // если ненулевые элементы под диагональю есть:
  // квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого
  // элемента
  if (std::sqrt(std::abs(s - A(j, i) * A(j, i))) >
      std::numeric_limits<ld>::min()) {
    // выбор знака слагаемого beta = sign(-x1)
    if (A(j, i) < 0)
      beta = std::sqrt(s);
    else
      beta = -std::sqrt(s);
    // вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
    mu = 1.0 / (beta * (beta - A(j, i)));
    // формируем вектор p
    for (int I = 0; I < A.getRows(); I++) {
      p(I, 0) = 0;
      if (I >= j) p(I, 0) = A(I, i);
    }
    // изменяем элемент, с которого начнётся обнуление
    p(j, 0) -= beta;
    // вычисляем новые компоненты матрицы A = ... * U2 * U1 * A
    for (int m = 0; m < A.getCols(); m++) {
      // произведение S = St * p
      s = 0;
      for (int n = j; n < A.getRows(); n++) {
        s += A(n, m) * p(n, 0);
      }
      s *= mu;
      // S = S - 2 * p * (St * p)^t / ||p||^2
      for (int n = j; n < A.getRows(); n++) {
        A(n, m) -= s * p(n, 0);
      }
    }
    // вычисляем новые компоненты матрицы U = ... * H2 * H1 * U
    for (int m = 0; m < A.getRows(); m++) {
      // произведение S = Ut * p
      s = 0;
      for (int n = j; n < A.getRows(); n++) {
        s += U(m, n) * p(n, 0);
      }
      s *= mu;
      // U = U - 2 * p * (Ut * p)^t / ||p||^2
      for (int n = j; n < A.getRows(); n++) {
        U(m, n) -= s * p(n, 0);
      }
    }
  }
}

void Householder_Row_Transform(Matrix& A, Matrix& V, int i, int j) {
  // вектор отражения
  Matrix p(A.getCols(), 1);
  // вспомогательные переменные
  ld s, beta, mu;
  // находим квадрат нормы строки для обнуления
  s = 0;
  for (int I = j; I < A.getCols(); I++) s += std::pow(A(i, I), 2);
  // если ненулевые элементы под диагональю есть:
  // квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого
  // элемента
  if (std::sqrt(std::abs(s - A(i, j) * A(i, j))) >
      std::numeric_limits<ld>::min()) {
    // выбор знака слагаемого beta = sign(-x1)
    if (A(i, j) < 0)
      beta = std::sqrt(s);
    else
      beta = -std::sqrt(s);
    // вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
    mu = 1.0 / (beta * (beta - A(i, j)));
    // формируем вектор p
    for (int I = 0; I < A.getCols(); I++) {
      p(I, 0) = 0;
      if (I >= j) p(I, 0) = A(i, I);
    }
    // изменяем диагональный элемент
    p(j, 0) -= beta;
    // вычисляем новые компоненты матрицы A = A * H1 * H2 ...
    for (int m = 0; m < A.getRows(); m++) {
      // произведение A * p
      s = 0;
      for (int n = j; n < A.getCols(); n++) {
        s += A(m, n) * p(n, 0);
      }
      s *= mu;
      // A = A - p * (A * p)^t
      for (int n = j; n < A.getCols(); n++) {
        A(m, n) -= s * p(n, 0);
      }
    }
    // вычисляем новые компоненты матрицы V = V * H1 * H2 * ...
    for (int m = 0; m < A.getCols(); m++) {
      // произведение V * p
      s = 0;
      for (int n = j; n < A.getCols(); n++) {
        s += V(m, n) * p(n, 0);
      }
      s *= mu;
      // V = V - p * (V * p)^t
      for (int n = j; n < A.getCols(); n++) {
        V(m, n) -= s * p(n, 0);
      }
    }
  }
}

void Givens_Delete_Elem_Down_Triangle(Matrix& A, Matrix& U, int I, int J) {
  ld help1, help2;

  // косинус, синус
  ld c = 0, s = 0;

  // если элемент не нулевой, то требуется поворот вектора
  if (std::abs(A(I, J)) > std::numeric_limits<ld>::min()) {
    help1 = std::sqrt(std::pow(A(I, J), 2) + std::pow(A(J, J), 2));
    c = A(J, J) / help1;
    s = A(I, J) / help1;

    // A_new = Gt * A
    for (int k = 0; k < A.getCols(); k++) {
      help1 = c * A(J, k) + s * A(I, k);
      help2 = c * A(I, k) - s * A(J, k);
      A(J, k) = help1;
      A(I, k) = help2;
    }
    // умножаем матрицу U на матрицу преобразования G справа: D = Qt * A * Q =>
    // Qt транспонируется для матрицы U
    for (int k = 0; k < U.getRows(); k++) {
      help1 = c * U(k, J) + s * U(k, I);
      help2 = c * U(k, I) - s * U(k, J);
      U(k, J) = help1;
      U(k, I) = help2;
    }
  }
  A(I, J) = 0;
}

void Givens_Delete_Elem_Up_Triangle(Matrix& A, Matrix& V, int I, int J) {
  ld help1, help2;

  // косинус, синус
  ld c = 0, s = 0;

  // если элемент не нулевой, то требуется поворот вектора
  if (std::abs(A(I, J)) > std::numeric_limits<ld>::min()) {
    help1 = std::sqrt(std::pow(A(I, J), 2) + std::pow(A(I, I), 2));
    c = A(I, I) / help1;
    s = -A(I, J) / help1;

    // A_new = A * Gt
    for (int k = 0; k < A.getRows(); k++) {
      help1 = c * A(k, I) - s * A(k, J);
      help2 = c * A(k, J) + s * A(k, I);
      A(k, I) = help1;
      A(k, J) = help2;
    }
    // умножаем матрицу V на матрицу преобразования Gt справа
    for (int k = 0; k < V.getRows(); k++) {
      help1 = c * V(k, I) - s * V(k, J);
      help2 = c * V(k, J) + s * V(k, I);
      V(k, I) = help1;
      V(k, J) = help2;
    }
  }
}

void Matrix::SVD(Matrix& U, Matrix& Sigma, Matrix& V, ld Reduction) const {
  // наименьшее измерение
  int Min_Size = std::min(rows_, cols_);
  // размеры нижней и верхней внешних диагоналей
  int Up_Size = Min_Size - 1, Down_Size = Min_Size - 1;
  // инициализация матрицы левых сингулярных векторов
  U.setRow(rows_);
  U.setCol(rows_);
  // матрица сингулярных чисел
  Sigma.setRow(rows_);
  Sigma.setCol(cols_);
  // инициализация матрицы правых сингулярных векторов
  V.setRow(cols_);
  V.setCol(cols_);
  // инициализация матриц для SVD
  for (int i = 0; i < rows_; i++) {
    U(i, i) = 1.0;
    for (int j = 0; j < cols_; j++) Sigma(i, j) = (*this)(i, j);
  }
  for (int i = 0; i < cols_; i++) V(i, i) = 1.0;

  //**************** этап I: бидиагонализация *************************

  for (int i = 0; i < Min_Size - 1; i++) {
    Householder_Col_Transform(Sigma, U, i, i);
    Householder_Row_Transform(Sigma, V, i, i + 1);
  }

  // ситуация M > N - строк больше => дополнительное умножение слева
  if (rows_ > cols_) {
    Householder_Col_Transform(Sigma, U, cols_ - 1, cols_ - 1);
    // нижняя побочная диагональ длиннее на 1
    Down_Size += 1;
  }

  // ситуация M < N - столбцов больше => дополнительное умножение справа
  if (rows_ < cols_) {
    Householder_Row_Transform(Sigma, V, rows_ - 1, rows_);
    // верхняя побочная диагональ длиннее на 1
    Up_Size += 1;
  }
  //**************** этап II: преследование ************
  //********* приведение к диагональному виду **********

  // для хранение изменяющихся элементов верхней диагонали
  ld Up[Up_Size];
  ld Down[Down_Size];
  // Initialize Up array with current upper diagonal elements
  for (int i = 0; i < Up_Size; i++) Up[i] = 0;
  // Initialize Down array with current lower diagonal elements
  for (int i = 0; i < Down_Size; i++) Down[i] = 0;
  // число неизменившихся элементов над главной диагональю
  int CountUpElements;
  // процедура преследования
  do {
    CountUpElements = 0;

    // обнуление верхней диагонали
    for (int i = 0; i < Up_Size; i++) {
      if (std::abs(Up[i] - Sigma(i, i + 1)) > std::numeric_limits<ld>::min()) {
        Up[i] = Sigma(i, i + 1);
        // Givens_Delete_Elem_Up_Triangle(Sigma, V, i, i + 1);
        Householder_Row_Transform(Sigma, V, i, i);
      } else
        CountUpElements++;
    }

    // обнуление нижней диагонали
    for (int i = 0; i < Down_Size; i++) {
      if (std::abs(Down[i] - Sigma(i + 1, i)) >
          std::numeric_limits<ld>::min()) {
        Down[i] = Sigma(i + 1, i);
        // Givens_Delete_Elem_Down_Triangle(Sigma, U, i + 1, i);
        Householder_Col_Transform(Sigma, U, i, i);
      }
    }
  } while (CountUpElements != Up_Size);
  //----------------------------------------
  // убираем отрицательные сингулярные числа
  // наименьшее измерение
  Min_Size = std::min(Sigma.getRows(), Sigma.getCols());
  // проверка сингулярных чисел на положительность
  for (int i = 0; i < Min_Size; i++) {
    if (Sigma(i, i) < 0) {
      Sigma(i, i) = -Sigma(i, i);
      for (int j = 0; j < U.getRows(); j++) U(j, i) = -U(j, i);
    }
  }
  //-----------------------------------------
  // сортируем по возрастанию сингулярные числа
  // сортировка сингулярных чисел
  for (int I = 0; I < Min_Size; I++) {
    ld Max_Elem = Sigma(I, I);
    int Index = I;
    for (int i = I + 1; i < Min_Size; i++) {
      if (Sigma(i, i) > Max_Elem) {
        Max_Elem = Sigma(i, i);
        Index = i;
      }
    }
    // найден наибольший элемент
    if (I != Index) {
      Sigma(Index, Index) = Sigma(I, I);
      Sigma(I, I) = Max_Elem;
      for (int j = 0; j < U.getRows(); j++) std::swap(U(j, I), U(j, Index));
      for (int j = 0; j < V.getRows(); j++) std::swap(V(j, I), V(j, Index));
    }
  }
  //-----------------------------------------
  // Reduction = Sigma(0, 0) * std::numeric_limits<ld>::epsilon();
  //  проверка на возможность редукции по сингулярным числам
  for (int i = 0; i < Min_Size; i++) {
    if (std::abs(Sigma(i, i)) < Reduction) {
      Min_Size = i;
      break;
    }
  }
  // редукция размерности матриц
  Sigma.setRow(Min_Size);
  Sigma.setCol(Min_Size);
  U.setCol(Min_Size);
  V.setCol(Min_Size);
}