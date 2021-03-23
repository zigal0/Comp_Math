#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

double FunctionG(double x);

double FunctionP(double x);

double FunctionF(double x);

double FunctionU(double x);

std::vector<double> ColumnF(int iterations, double step, double icLeft, double icRight);

std::vector<std::vector<double>> FillMatrix(int iterations, double step);

int main() {
    // introduction
    std::cout << "A solution to a boundary value problem: - (g(x) * u'(x))' + p(x) * u(x) = f(x)" << std::endl;
    std::cout << "with initial conditions: u'(leftBoundary) = a & u'(rightBoundary) = b" << std::endl;
    // set initial data
    double leftBoundary = 0;
    double rightBoundary = 1;
    int iterations = 20;
    double step = (rightBoundary - leftBoundary) / iterations;

    //// initial conditions for tests
    //// test#1
//    double icLeft = exp(0);
//    double icRight = exp(1);
    //// test#2
//    double icLeft = cos(leftBoundary);
//    double icRight = cos(rightBoundary);
    //// test#3
//    double icLeft = 1 / (cos(leftBoundary) * cos(leftBoundary));
//    double icRight = 1 / (cos(rightBoundary) * cos(rightBoundary));
    //// test#4
    double icLeft = exp(leftBoundary) + 1;
    double icRight = exp(rightBoundary) + 1;

    // creation different structures for task
    std::vector<double> F = ColumnF(iterations, step, icLeft, icRight);
    std::vector<std::vector<double>> M = FillMatrix(iterations, step);
    std::vector<double> U_APPROX(iterations + 1);

    // check diagonally dominant
    bool check = true;
    for (int i = 0; i < iterations + 1; i++) {
        if (std::abs(M[i][1]) < std::abs(M[i][0]) + std::abs(M[i][2])) {
            std::cout << "Not possible" << std::endl;
            check = false;
            break;
        }
    }
    if (check) {
        // straight
        for (int i = 1; i < iterations + 1; i++) {
            M[i][1] = M[i][1] - (M[i][0] * M[i - 1][2]) / M[i - 1][1];
            F[i] = F[i] - (M[i][0] * F[i - 1]) / M[i - 1][1];
            M[i][0] = 0;
        }
        // reverse
        U_APPROX[iterations] = F[iterations] / M[iterations][1];
        for (int i = iterations - 1; i >= 0; i--) {
            U_APPROX[i] = (F[i] - M[i][2] * U_APPROX[i + 1]) / M[i][1];
        }
        // cubic norm + real solution
        double cubNorm = 0;
        std::vector<double> U_REAL(iterations + 1);
        for (int j = 0; j < iterations + 1; j++) {
            U_REAL[j] = FunctionU(j * step);
            if (std::abs(U_REAL[j] - U_APPROX[j]) > cubNorm) {
                cubNorm = std::abs(U_REAL[j] - U_APPROX[j]);
            }
        }

        // Results:
        std::cout << "x_cur" << std::setw(21) << "U_REAl" << std::setw(22) << "U_APPROX" << std::endl;
        for (int j = 0; j < iterations + 1; j++) {
            std::cout << std::fixed << std::setprecision(10)
                      << step * j
                      << std::setw(20) << U_REAL[j]
                      << std::setw(20) << U_APPROX[j]
                      << std::endl;
        }
        std::cout << "Cubic norm = " << cubNorm << std::endl;
    }
    return 0;
}

std::vector<double> ColumnF(int iterations, double step, double icLeft, double icRight) {
    std::vector<double> res(iterations + 1);
    // first value
    res[0] = FunctionF(0) - FunctionG(-(step / 2)) * (icLeft * 2) / step;
    // last value
    res[iterations] = FunctionF(iterations * step) + FunctionG(iterations * step + step / 2) * (icRight * 2) / step;
    // other values
    for (int i = 1; i < iterations; i++) {
        res[i] = FunctionF(i * step);
    }
    return res;
}

std::vector<std::vector<double>> FillMatrix(int iterations, double step) {
    std::vector<std::vector<double>> M(iterations + 1, std::vector<double>(3));
    // 1st row
    M[0][0] = 0;
    M[0][1] = FunctionP(0) + (FunctionG(-(step / 2)) + FunctionG(step / 2)) / (step * step);
    M[0][2] = -(FunctionG(-(step / 2)) + FunctionG(step / 2)) / (step * step);
    // last row
    M[iterations][0] =
            -(FunctionG(step * iterations + step / 2) + FunctionG(step * iterations - step / 2)) / (step * step);
    M[iterations][1] = FunctionP(iterations * step) +
                       (FunctionG(step * iterations + step / 2) + FunctionG(step * iterations - step / 2)) /
                       (step * step);
    M[iterations][2] = 0;
    // other rows
    for (int i = 1; i < iterations; i++) {
        M[i][0] = -FunctionG(step * i - step / 2) / (step * step);
        M[i][1] =
                FunctionP(i * step) + (FunctionG(step * i - step / 2) + FunctionG(step * i + step / 2)) / (step * step);
        M[i][2] = -FunctionG(step * i + step / 2) / (step * step);
    }
    return M;
}

//// TEST

////// test#1 (e^x)
//double FunctionG(double x) {
//    return (x * x + 1);
//}
//
//double FunctionP(double x) {
//    return (x + 1);
//}
//
//double FunctionF(double x) {
//    return -x * exp(x) * (x + 1);
//}
//
//double FunctionU(double x) {
//    return exp(x);
//}

//// test#2 (sin(x) + 1)
//double FunctionG(double x) {
//    return (x * x + 1);
//}
//
//double FunctionP(double x) {
//    return (x + 1);
//}
//
//double FunctionF(double x) {
//    return (-2 * x * cos(x) + (x * x + 1) * sin(x) + (x + 1) * (sin(x) + 1));
//}
//
//double FunctionU(double x) {
//    return (sin(x) + 1);
//}

//// test#3 (tan(x) + 1)
//double FunctionG(double x) {
//    return exp(x);
//}
//
//double FunctionP(double x) {
//    return exp(x);
//}
//
//double FunctionF(double x) {
//    return (exp(x) * ((-1 - 2 * tan(x)) / (cos(x) * cos(x)) + tan(x) + 1));
//}
//
//double FunctionU(double x) {
//    return (tan(x) + 1);
//}

//// test#4 (tan(x) + 1)
double FunctionG(double x) {
    return x * x + 1;
}

double FunctionP(double x) {
    return x + 1;
}

double FunctionF(double x) {
    return (-2 * x * (exp(x) + 1) - (x * x + 1) * exp(x) + (x+1) * (exp(x) + x));
}

double FunctionU(double x) {
    return (exp(x) + x);
}