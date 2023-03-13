#include "Includes.h"

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