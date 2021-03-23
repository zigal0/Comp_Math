#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

struct Parameters {
    double l = 1;
    double lambda = 0.1;
    double g = 9.8;
    double u_0 = 0.01;
    double v_0 = 0.05;
    int n_steps = 20;
    int numberOfPeriods = 3;
};

int main() {
    // introduction
    std::cout << "Modulation of vibration: l * u''(t) + 2 * lambda * u'(t) + g * sin(u(t)) = 0" << std::endl;
    std::cout << "with initial conditions: u(0) = u_0 & u'(0) = v_0" << std::endl;
    std::string cmd;
    Parameters param;
    do {
        std::cout << " - Enter 'RUN' if you want to run program;" << std::endl;
        std::cout << " - Enter 'SHOW' if you want to see parameters;" << std::endl;
        std::cout << " - Enter 'CHANGE' if you want to change parameters;" << std::endl;
        std::cin >> cmd;
        if (cmd == "CHANGE") {
            std::cout << "Enter next parameters: " << std::endl;
            // parameters
            std::cout << "length =";
            std::cin >> param.l;
            std::cout << "lambda =";
            std::cin >> param.lambda;
            std::cout << "g =";
            std::cin >> param.g;
            std::cout << "u_0 =";
            std::cin >> param.u_0;
            std::cout << "v_0 =";
            std::cin >> param.v_0;
            do {
                std::cout << "numberOfPeriods =";
                std::cin >> param.numberOfPeriods;
                if (param.numberOfPeriods <= 0) {
                    std::cout << "Wrong number of periods, try other" << std::endl;
                }
            } while (param.numberOfPeriods <= 0);
            do {
                std::cout << "n_steps =";
                std::cin >> param.n_steps;
                if (param.n_steps <= 0) {
                    std::cout << "Wrong quantity steps, try other" << std::endl;
                }
            } while (param.n_steps <= 0);
        } else if (cmd == "SHOW") {
            std::cout << "__________________________________________" << std::endl;
            std::cout << "length = " << param.l << std::endl;
            std::cout << "lambda = " << param.lambda << std::endl;
            std::cout << "g = " << param.g << std::endl;
            std::cout << "u_0 = " << param.u_0 << std::endl;
            std::cout << "v_0 = " << param.v_0 << std::endl;
            std::cout << "n_steps = " << param.n_steps << std::endl;
            std::cout << "numberOfPeriods = " << param.numberOfPeriods << std::endl;
            std::cout << "__________________________________________" << std::endl;
        } else if (cmd != "RUN") {
            std::cout << "Wrong command, try other" << std::endl;
        }
    } while (cmd != "RUN");

    // parameters for solution
    double k = param.lambda / param.l;
    double w_0 = sqrt(param.g / param.l);
    double T = (2 * M_PI) / w_0;
    double tay = T / param.n_steps;
    double w = w_0 * w_0 - k * k;
    int n = param.n_steps * param.numberOfPeriods;
    std::vector<std::vector<double>> res(4, std::vector<double>(n));
    for (int j = 0; j < n; j++) {
        res[0][j] = tay * j;
    }

    // AN_approach
    if (param.u_0 < 0.1) {
        double u_AN;
        if (w > 0) {
//            std::cout << "weak attenuation" << std::endl;
            w = sqrt(w);
            double a = param.u_0;
            double b = (param.v_0 + param.u_0 * k) / w;
            double t_cur = 0;
            for (int j = 0; j < n; j++) {
                u_AN = exp(-k * t_cur) * (a * cos(w * t_cur) + b * sin(w * t_cur));
                res[1][j] = u_AN;
                t_cur += tay;
            }
        } else if (w < 0) {
//            std::cout << "strong attenuation" << std::endl;
            w = sqrt(-w);
            double a = 0.5 * (param.u_0 + (param.v_0 + k * param.u_0) / w);
            double b = 0.5 * (param.u_0 - (param.v_0 + k * param.u_0) / w);
            double t_cur = 0;
            for (int j = 0; j < n; j++) {
                u_AN = exp(-k * t_cur) * (a * exp(w * t_cur) + b * exp(-w * t_cur));
                res[1][j] = u_AN;
                t_cur += tay;
            }
        } else {
//            std::cout << "degeneration" << std::endl;
            double a = param.u_0;
            double b = param.v_0 + param.u_0 * k;
            double t_cur = 0;
            for (int j = 0; j < n; j++) {
                u_AN = exp(-k * t_cur) * (a + b * t_cur);
                res[1][j] = u_AN;
                t_cur += tay;
            }
        }
    }

    // EE_approach
    std::vector<double> vel(n);
    res[2][0] = param.u_0;
    vel[0] = param.v_0;
    for (int j = 0; j < n - 1; j++) {
        res[2][j + 1] = res[2][j] + tay * vel[j];
        vel[j + 1] = vel[j] - tay * (2 * k * vel[j] + w_0 * w_0 * sin(res[2][j]));
    }

    // RK_approach
    res[3][0] = param.u_0;
    vel[0] = param.v_0;
    double k_coef[2][4] = {0};
    for (int j = 0; j < n - 1; j++) {
        k_coef[0][0] = vel[j];
        k_coef[1][0] = -2 * k * vel[j] - w_0 * w_0 * sin(res[3][j]);

        k_coef[0][1] = vel[j] + tay * k_coef[1][0] / 2;
        k_coef[1][1] = -2 * k * (vel[j] + tay * k_coef[1][0] / 2) - w_0 * w_0 * sin(res[3][j] + tay * k_coef[0][0] / 2);

        k_coef[0][2] = vel[j] + tay * k_coef[1][1] / 2;
        k_coef[1][2] = -2 * k * (vel[j] + tay * k_coef[1][1] / 2) - w_0 * w_0 * sin(res[3][j] + tay * k_coef[0][1] / 2);

        k_coef[0][3] = vel[j] + tay * k_coef[1][2];
        k_coef[1][3] = -2 * k * (vel[j] + tay * k_coef[1][2]) - w_0 * w_0 * sin(res[3][j] + tay * k_coef[0][2]);

        // Coordinate
        res[3][j + 1] = res[3][j] + tay * (k_coef[0][0] + 2 * (k_coef[0][1] + k_coef[0][2]) + k_coef[0][3]) / 6;

        // Velocity
        vel[j + 1] = vel[j] + tay * (k_coef[1][0] + 2 * (k_coef[1][1] + k_coef[1][2]) + k_coef[1][3]) / 6;
    }
    // Results:
    std::cout  << "t_cur" << std::setw(20) << "U_AN" << std::setw(20) << "U_EE" << std::setw(20) << "U_RK" << std::endl;
    for (int j = 0; j < n; j++) {
        std::cout << std::fixed << std::setprecision(10)
                  << res[0][j]
                  << std::setw(20) << res[1][j]
                  << std::setw(20) << res[2][j]
                  << std::setw(20) << res[3][j]
                  << std::endl;
    }
    return 0;
}
