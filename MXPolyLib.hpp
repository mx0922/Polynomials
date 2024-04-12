#pragma once
#ifndef MXPOLYLIB_HPP
#define MXPOLYLIB_HPP

#include <vector>
#include <cmath>
#include <algorithm>

// function: clamp the s in [0, 1]
inline double clamp01(double s){
    return std::max(0.0, std::min(1.0, s));
}

// MX_Poly5 class:
class MX_Poly5 {
public:
    MX_Poly5(double t, double p0, double v0, double a0, double p1, double v1, double a1)
        : T(t), P0(p0), V0(v0), A0(a0), P1(p1), V1(v1), A1(a1), S(6) {
        autoGen_fiveOrderPolyCoeff();
    }

    void getFiveOrderPoly(double s, double& pos, double& vel, double& acc) const {
        double t = clamp01(s) * T;

        pos = S[0] + S[1] * t + S[2] * t * t + S[3] * std::pow(t, 3) + S[4] * std::pow(t, 4) + S[5] * std::pow(t, 5);
        vel = S[1] + 2 * S[2] * t + 3 * S[3] * std::pow(t, 2) + 4 * S[4] * std::pow(t, 3) + 5 * S[5] * std::pow(t, 4);
        acc = 2 * S[2] + 6 * S[3] * t + 12 * S[4] * std::pow(t, 2) + 20 * S[5] * std::pow(t, 3);
    }

private:
    double T;
    double P0;
    double V0;
    double A0;
    double P1;
    double V1;
    double A1;
    std::vector<double> S;

    void autoGen_fiveOrderPolyCoeff() {
        S[0] = P0;
        S[1] = V0;
        S[2] = A0 / 2.0;

        double t2 = T * T;
        double t3 = A1 * t2;
        double t4 = A0 * t2 * 3.0;
        double t5 = -t3;

        S[3] = 1.0 / std::pow(T, 3) * (P0 * 2.0e+1 - P1 * 2.0e+1 + t4 + t5 + T * V0 * 1.2e+1 + T * V1 * 8.0) * (-1.0 / 2.0);
        S[4] = (1.0 / std::pow(t2, 2) * (P0 * 3.0e+1 - P1 * 3.0e+1 - t3 * 2.0 + t4 + T * V0 * 1.6e+1 + T * V1 * 1.4e+1)) / 2.0;
        S[5] = 1.0 / std::pow(T, 5) * (P0 * 1.2e+1 - P1 * 1.2e+1 + t5 + T * V0 * 6.0 + T * V1 * 6.0 + A0 * t2) * (-1.0 / 2.0);
    }
};

// MX_Poly6 class:
class MX_Poly6 {
public:
    MX_Poly6(double t, double p0, double v0, double a0, double pm, double p1, double v1, double a1)
        : T(t), P0(p0), V0(v0), A0(a0), Pm(pm), P1(p1), V1(v1), A1(a1), S(7) {
        autoGen_sixOrderPolyCoeff();
    }

    void getSixOrderPoly(double s, double& pos, double& vel, double& acc) const {
        double t = clamp01(s) * T;

        pos = S[0] + S[1] * t + S[2] * std::pow(t, 2) + S[3] * std::pow(t, 3) + S[4] * std::pow(t, 4) +
              S[5] * std::pow(t, 5) + S[6] * std::pow(t, 6);

        vel = S[1] + 2 * S[2] * t + 3 * S[3] * std::pow(t, 2) + 4 * S[4] * std::pow(t, 3) +
              5 * S[5] * std::pow(t, 4) + 6 * S[6] * std::pow(t, 5);

        acc = 2 * S[2] + 6 * S[3] * t + 12 * S[4] * std::pow(t, 2) + 20 * S[5] * std::pow(t, 3) +
              30 * S[6] * std::pow(t, 4);
    }

private:
    double T;
    double P0;
    double V0;
    double A0;
    double Pm;
    double P1;
    double V1;
    double A1;
    std::vector<double> S;

    void autoGen_sixOrderPolyCoeff() {
    S[0] = P0;
    S[1] = V0;
    S[2] = A0 / 2.0;

    double t2 = T * T;
    double t4 = Pm * 3.84e+2;
    double t3 = A1 * t2;

    S[3] = 1.0 / std::pow(T, 3) * (P0 * 8.4e+1 + P1 * 4.4e+1 - Pm * 1.28e+2 + t3 +
                                   T * V0 * 3.2e+1 - T * V1 * 1.2e+1 + A0 * t2 * 5.0) * (-1.0 / 2.0);

    S[4] = (1.0 / std::pow(t2, 2) * (P0 * 2.22e+2 + P1 * 1.62e+2 + t3 * 4.0 - t4 +  // Fix: Corrected the sign here
                                     T * V0 * 7.6e+1 - T * V1 * 4.6e+1 + A0 * t2 * 9.0)) / 2.0;

    S[5] = 1.0 / std::pow(T, 5) * (P0 * 2.04e+2 + P1 * 1.8e+2 + t3 * 5.0 - t4 +
                                   T * V0 * 6.6e+1 - T * V1 * 5.4e+1 + A0 * t2 * 7.0) * (-1.0 / 2.0);

    S[6] = 1.0 / std::pow(t2, 3) * (P0 * 3.2e+1 + P1 * 3.2e+1 - Pm * 6.4e+1 + t3 +
                                     T * V0 * 1.0e+1 - T * V1 * 1.0e+1 + A0 * t2);
    }
};

#endif