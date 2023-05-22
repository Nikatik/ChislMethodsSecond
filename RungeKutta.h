#include "Includes.h"

// Runge-Kutta in general form (one step)
void RK (
    double, double*, double*, double*, double*, double*, double, unsigned int,
    double**, double**,
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double, double (*) (double, double, double, double, double, double), bool);

// The chord method
void diagonal (
    double, double, double, double, double, double, double*, double*, double*,
    double*, double*, double*, unsigned int, double**, double**, double,
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double, double (*) (double, double, double, double, double, double),
    double (*) (double, double, double, double, double, double), bool);

// Adaptive Runge-Kutta
double astep (
    double, double*, double*, double*, double*, double*, long long unsigned*,
    long long unsigned*, unsigned int, unsigned int, double**, double**, double,
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double (*) (double, double, double, double, double, double,
                double (*) (double, double, double, double, double, double)),
    double, bool, bool);

void RKq (
    __float128, __float128*, __float128*, __float128*, __float128*, __float128*, __float128, unsigned int,
    __float128**, __float128**,
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128, __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128), bool);

// The chord method
void diagonalq (
    __float128, __float128, __float128, __float128, __float128, __float128, __float128*, __float128*, __float128*,
    __float128*, __float128*, __float128*, unsigned int, __float128**, __float128**, __float128,
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128, __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128), bool);

// Adaptive Runge-Kutta
__float128 astepq (
    __float128, __float128*, __float128*, __float128*, __float128*, __float128*, long long unsigned*,
    long long unsigned*, unsigned int, unsigned int, __float128**, __float128**, __float128,
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128,
                __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128, bool, bool);
