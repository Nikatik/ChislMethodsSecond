#include <math.h>
#include <stdio.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

static double
    x_d (double t, double x, double y, double p_x, double p_y, double mult,
         double control (double, double, double, double, double, double))
{
    return y;
}

static double
    y_d (double t, double x, double y, double p_x, double p_y, double mult,
         double control (double, double, double, double, double, double))
{
    return control (t, x, y, p_x, p_y, mult);
}

static double
    px_d (double t, double x, double y, double p_x, double p_y, double mult,
          double control (double, double, double, double, double, double))
{
    return (-mult / 2.) * sin (mult * x) *
           pow (control (t, x, y, p_x, p_y, mult), 2);
}

static double
    py_d (double t, double x, double y, double p_x, double p_y, double mult,
          double control (double, double, double, double, double, double))
{
    return -p_x - y;
}

static double
    B_d (double t, double x, double y, double p_x, double p_y, double mult,
         double contol (double, double, double, double, double, double))
{
    return pow (x, 2.) + pow (p_y, 2.) / 4.;
}

static double
    test_x_d (double t, double x, double y, double p_x, double p_y, double mult,
              double contol (double, double, double, double, double, double))
{
    return -y;
}

static double
    test_y_d (double t, double x, double y, double p_x, double p_y, double mult,
              double contol (double, double, double, double, double, double))
{
    return x;
}

static double test_px_d (double t, double x, double y, double p_x, double p_y,
                         double mult,
                         double contol (double, double, double, double, double,
                                        double))
{
    return -p_y;
}

static double test_py_d (double t, double x, double y, double p_x, double p_y,
                         double mult,
                         double contol (double, double, double, double, double,
                                        double))
{
    return p_x;
}

static double
    test_B_d (double t, double x, double y, double p_x, double p_y, double mult,
              double contol (double, double, double, double, double, double))
{
    return sqrt (pow (x - p_x, 2) + pow (y - p_y, 2));
}

#pragma GCC diagnostic pop
