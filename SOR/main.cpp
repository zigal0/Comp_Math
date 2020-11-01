#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>

class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    int num_rows;
    int num_cols;
public:
    // Default constructor
    Matrix() {
        num_rows = 0;
        num_cols = 0;
    }

    // Constructor
    Matrix(int rows, int cols) {
        num_rows = 0;
        num_cols = 0;
        Reset(rows, cols);
    }

    // Change number of rows and number of columns
    void Reset(int rows, int cols) {
        if (rows < 0) {
            throw std::out_of_range("Num_rows must be >= 0");
        }
        if (cols < 0) {
            throw std::out_of_range("Num_cols must be >= 0");
        }
        if (cols == 0 || rows == 0) {
            cols = rows = 0;
        }
        num_rows = rows;
        num_cols = cols;
        matrix.assign(num_rows, std::vector<double>(num_cols));
    }

    // Return element with (row, col)
    [[nodiscard]] double At(const int &row, const int &col) const {
        return matrix.at(row).at(col);
    }

    // Return reference of element with (row, col)
    [[nodiscard]] double &At(const int &row, const int &col) {
        return matrix.at(row).at(col);
    }

    // Return number of rows
    [[nodiscard]] int GetNumRows() const {
        return num_rows;
    }

    // Return number of cols
    [[nodiscard]] int GetNumColumns() const {
        return num_cols;
    }
};

//// Normal input
//std::istream &operator>>(std::istream &in, Matrix &m) {
//    int cols = 0, rows = 0;
//    in >> rows >> cols;
//    m.Reset(rows, cols);
//    for (int row = 0; row < rows; row++) {
//        for (int col = 0; col < cols; col++) {
//            in >> m.At(row, col);
//        }
//    }
//    return in;
//}

// Input for task
std::istream &operator>>(std::istream &in, Matrix &m) {
    int size = 0;
    in >> size;
    m.Reset(size, size);
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            in >> m.At(row, col);
        }
    }
    return in;
}

// Output
std::ostream &operator<<(std::ostream &out, const Matrix &m) {

    for (int row = 0; row < m.GetNumRows(); row++) {
        for (int col = 0; col < m.GetNumColumns(); col++) {
            if (col > 0) {
                out << ' ';
            }
//            out << std::fixed << std::setprecision(2) << m.At(row, col);
            out << m.At(row, col);
        }
        out << std::endl;
    }
    return out;
}

//// Additional function for Matrix
//// Check for equality
//bool operator==(const Matrix &m1, const Matrix &m2) {
//    if (m1.GetNumRows() != m2.GetNumRows() || m1.GetNumColumns() != m2.GetNumColumns()) {
//        return false;
//    }
//    for (int row = 0; row < m1.GetNumRows(); row++) {
//        for (int col = 0; col < m1.GetNumColumns(); col++) {
//            if (m1.At(row, col) != m2.At(row, col)) {
//                return false;
//            }
//        }
//    }
//    return true;
//}
//
//// Multiplication matrix by number №1
//Matrix operator*(double factor, const Matrix &m) {
//    if (m.GetNumColumns() == 0) {
//        return m;
//    }
//    Matrix res(m.GetNumRows(), m.GetNumColumns());
//    for (int row = 0; row < m.GetNumRows(); row++) {
//        for (int col = 0; col < m.GetNumColumns(); col++) {
//            res.At(row, col) = factor * m.At(row, col);
//        }
//    }
//    return res;
//}
//
//// Multiplication matrix by number №2
//Matrix operator*(const Matrix &m, double factor) {
//    return factor * m;
//}
//
//// Division of matrix by number
//Matrix operator/(const Matrix &m, double factor) {
//    if (factor == 0) {
//        throw std::invalid_argument("Division by zero");
//    }
//    return m * (1 / factor);
//}
//
//// Sum of two matrices
//Matrix operator+(const Matrix &m1, const Matrix &m2) {
//    if (m1.GetNumRows() != m2.GetNumRows()) {
//        throw std::invalid_argument("Mismatched number of rows");
//    }
//
//    if (m1.GetNumColumns() != m2.GetNumColumns()) {
//        throw std::invalid_argument("Mismatched number of columns");
//    }
//    Matrix res(m1.GetNumRows(), m1.GetNumColumns());
//    for (int row = 0; row < m1.GetNumRows(); row++) {
//        for (int col = 0; col < m1.GetNumColumns(); col++) {
//            res.At(row, col) = m1.At(row, col) + m2.At(row, col);
//        }
//    }
//    return res;
//}
//
//// Difference of two matrices
//Matrix operator-(const Matrix &m1, const Matrix &m2) {
//    Matrix m2_ = -1 * m2;
//    return m1 + m2_;
//}
//
//// Multiplication of two matrices
//Matrix operator*(const Matrix &m1, const Matrix &m2) {
//    if (m1.GetNumColumns() != m2.GetNumRows()) {
//        throw std::invalid_argument("Matrix multiplication is impossible");
//    }
//    if (m1.GetNumColumns() == 0 || m2.GetNumColumns() == 0) {
//        throw std::invalid_argument("Matrix multiplication is impossible");
//    }
//
//    Matrix res(m1.GetNumRows(), m2.GetNumColumns());
//    int i = 0, j = 0;
//    for (int row = 0; row < m1.GetNumRows(); row++) {
//        for (int lol = 0; lol < m2.GetNumColumns(); lol++) {
//            res.At(i, j) = 0;
//            for (int col = 0; col < m1.GetNumColumns(); col++) {
//                res.At(i, j) += m1.At(row, col) * m2.At(col, lol);
//            }
//            j++;
//        }
//        j = 0;
//        i++;
//    }
//    return res;
//}
//
//// Matrix Transpose
//Matrix Transpose(const Matrix &m) {
//    Matrix m_t(m.GetNumColumns(), m.GetNumRows());
//    for (int col = 0; col < m.GetNumColumns(); col++) {
//        for (int row = 0; row < m.GetNumRows(); row++) {
//            m_t.At(col, row) = m.At(row, col);
//        }
//    }
//    return m_t;
//}

// Multiplication vector by matrix №1
std::vector<double> operator*(const Matrix &m, std::vector<double> &v) {
    if (v.size() != m.GetNumColumns() || v.size() != m.GetNumRows()) {
        throw std::out_of_range("Invalid operation!");
    }
    std::vector<double> res(v.size(), 0);
    if (m.GetNumColumns() == 0) {
        return res;
    }
    for (int row = 0; row < m.GetNumRows(); row++) {
        for (int col = 0; col < m.GetNumColumns(); col++) {
            res.at(row) += m.At(row, col) * v.at(col);
        }
    }
    return res;
}

// Multiplication vector by matrix №2
std::vector<double> operator*(std::vector<double> &v, const Matrix &m) {
    return m * v;
}

// Output for vector
std::ostream &operator<<(std::ostream &out, const std::vector<double> &v) {
    for (double i : v) {
        std::cout << i << ' ';
    }
    return out;
}

//void next_iteration_z(double *x_pr, double *x_next, int n, double *matrix, double *d, double om) {
//    int i, j;
//    double s1, s2;
//    for (i = 0; i < n; i++) {
//        s1 = 0;
//        s2 = 0;
//        for (j = 0; j < i; j++) {
//            c[n * i + j] = -matrix[n * i + j] * om / matrix[n * i + i];
//            s1 = s1 + c[n * i + j] * x_next[j];
//        }
//        for (j = i + 1; j < n; j++) {
//            c[n * i + j] = -matrix[n * i + j] * om / matrix[n * i + i];
//            s2 = s2 + c[n * i + j] * x_pr[j];
//        }
//        d[i] = f[i] * om / matrix[n * i + i];
//        x_next[i] = s1 + s2 + d[i] - x_pr[i] * (om - 1);
//    }
//}

//// Finding the lower triangular matrix
//void Lower_triangular(const Matrix &A, Matrix &L) {
//
//    if (A.GetNumRows() != A.GetNumColumns()) {
//        throw std::invalid_argument("Matrix is not square!");
//    }
//    int size = A.GetNumColumns();
//    if (A.At(0, 0) <= 0) {
//        throw std::invalid_argument("Impossible sqrt!");
//    } else {
//        L.At(0, 0) = sqrt(A.At(0, 0));
//    }
//
//    for (int row = 1; row < size; row++) {
//        L.At(row, 0) = A.At(row, 0) / L.At(0, 0);
//    }
//
//    for (int i = 1; i < size; i++) {
//
//        double res = 0;
//        for (int p = 0; p <= i - 1; p++) {
//            res += L.At(i, p) * L.At(i, p);
//        }
//        res = A.At(i, i) - res;
//        if (res <= 0) {
//            throw std::invalid_argument("Impossible sqrt!");
//        } else {
//            L.At(i, i) = sqrt(res);
//            res = 0;
//        }
//        for (int j = i + 1; j < size; j++) {
//            for (int p = 0; p <= i - 1; p++) {
//                res += L.At(i, p) * L.At(j, p);
//            }
//            L.At(j, i) = (A.At(j, i) - res) / L.At(i, i);
//            res = 0;
//        }
//    }
//}
//
//// Finding y
//void Finding_y(const Matrix &L, std::vector<double> &y, std::vector<double> &b) {
//    int n = L.GetNumRows();
//    for (int i = 0; i < n; i++) {
//
//        double res = 0;
//        for (int j = 0; j <= i - 1; j++) {
//            res += y[j] * L.At(i, j);
//        }
//
//        if (L.At(i, i) == 0 && b[i] - res != 0) {
//            throw std::invalid_argument("Impossible division!");
//        } else {
//            y[i] = (b[i] - res) / L.At(i, i);
//        }
//    }
//}
//
//// Finding x_alg
//void Finding_x(const Matrix &L, std::vector<double> &x, std::vector<double> &y) {
//    int n = y.size();
//    for (int i = n - 1; i >= 0; i--) {
//
//        double res = 0;
//        for (int j = i + 1; j < n; j++) {
//            res = res + x[j] * L.At(j, i);
//        }
//
//        if (L.At(i, i) == 0 && y[i] - res != 0) {
//            throw std::invalid_argument("Impossible division!");
//        } else {
//            x[i] = (y[i] - res) / L.At(i, i);
//        }
//    }
//}
//
//// Finding Norm of Difference
//double Norm_of_difference(const std::vector<double> &x, const std::vector<double> &y) {
//    std::vector<double> diff(x.size());
//    int size = x.size();
//    double norm = 0;
//    for (int i = 0; i < size; i++) {
//        diff[i] = y[i] - x[i];
//        norm += diff[i] * diff[i];
//    }
//    norm = sqrt(norm);
//    return norm;
//}
//
//// Gauss–Seidel method
//void Next_iteration_GS(std::vector<double> &x, const Matrix &A, std::vector<double> &f) {
//    int i, j;
//    int n = A.GetNumRows();
//    double s1, s2;
//    for (i = 0; i < n; i++) {
//        s1 = 0;
//        s2 = 0;
//        for (j = 0; j < i; j++) {
//            s1 = s1 - (A.At(i, j) / A.At(i, i)) * x[j];
//        }
//        for (j = i + 1; j < n; j++) {
//            s2 = s2 - (A.At(i, j) / A.At(i, i)) * x[j];
//        }
//        x[i] = s1 + s2 + f[i] / A.At(i, i);
//    }
//}

//Successive over-relaxation method
void Next_iteration_SOR(const Matrix &A, std::vector<double> &x, std::vector<double> &f) {
    int k, j;
    double sum1, sum2;
    int n = A.GetNumRows();
    for (k = 0; k < n; k++) {
        sum1 = 0;
        sum2 = 0;
        for (j = 0; j < k; j++) {
            sum1 = sum1 - A.At(k, j) / A.At(k, k) * x[j];
        }
        for (j = k + 1; j < n; j++) {
            sum2 = sum2 - A.At(k, j) / A.At(k, k) * x[j];
        }
        x[k] = (sum1 + sum2 + f[k] / A.At(k, k)) * 1.5  - 0.5 * x[k];
    }
}

// Residual
double Residual( Matrix& A, std::vector<double> &x, const std::vector<double> &f) {
    std::vector<double> diff = A * x;
    int size = x.size();
    double norm = 0;
    for (int i = 0; i < size; i++) {
        diff[i] = f[i] - diff[i];
        norm += diff[i] * diff[i];
    }
    norm = sqrt(norm);
    return norm;
}


int main() {
    try {
        Matrix A;
        std::ifstream fin("Tests/matrix_2.txt");
        if (!fin.is_open())
            std::cout << "File cannot be opened!\n";
        else {
            fin >> A;
        }
        //Zero approximation
        std::vector<double> x{0, 0, 0, 0};
        // Exact solution
        std::vector<double> exsol {0, 0, 0, 0};
        srand(time(NULL));
        for (auto &i : exsol) {
            i = -100 + rand() % 200;
        }

        std::vector<double> f = exsol * A;

        // Results
        long int i = 1;
        double res = 0;
        do {
            std::cout << "Step № " << i << ':' << std::endl;
            i++;
            Next_iteration_SOR(A, x, f);
            std::cout << "Iteration : " << x << std::endl;
            res = Residual(A, x, f);
            std::cout << "Residual : " << res << std::endl;
        } while (res > 1e-4);
        std::cout << "Exact solution : " << exsol << std::endl;
    } catch (std::exception &ex) {
        std::cout << ex.what() << std::endl;
    }
    return 0;
}
