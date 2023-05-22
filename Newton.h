#include "Includes.h"

// Discrepancy counting
__float128 errorq (
    __float128, __float128, __float128, __float128, __float128, __float128, __float128*, unsigned int,
    unsigned int, __float128**, __float128**, __float128,
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
    __float128, bool);

// The shooting method
int shooting_methodq (
    unsigned int, __float128, __float128, __float128, __float128*, unsigned int, unsigned int,
    __float128**, __float128**, __float128,
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
    __float128);

// Discrepancy counting
double error (
    double, double, double, double, double, double, double*, unsigned int,
    unsigned int, double**, double**, double,
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
    double, bool);

// The shooting method
int shooting_method (
    unsigned int, double, double, double, double*, unsigned int, unsigned int,
    double**, double**, double,
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
    double);
