#include <math.h>
#include <quadmath.h>
#include <stdio.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

static __float128
    x_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
         __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return y;
}

static __float128
    y_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
         __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return control (t, x, y, p_x, p_y, mult);
}

static __float128
    px_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
          __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return (-mult / 2) * sinq (mult * x) *
           powq (control (t, x, y, p_x, p_y, mult), 2);
}

static __float128
    py_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
          __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return -p_x - y;
}

static __float128
    B_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
         __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return cosq (mult * x) * powq (control (t, x, y, p_x, p_y, mult), 2) -
           powq (y, 2);
}

static __float128
    test_x_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
              __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return -y;
}

static __float128
    test_y_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
              __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return x;
}

static __float128 test_px_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
                         __float128 mult,
                         __float128 control (__float128, __float128, __float128, __float128, __float128,
                                         __float128))
{
    return -p_y;
}

static __float128 test_py_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
                         __float128 mult,
                         __float128 control (__float128, __float128, __float128, __float128, __float128,
                                         __float128))
{
    return p_x;
}

static __float128
    test_B_dq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y, __float128 mult,
              __float128 control (__float128, __float128, __float128, __float128, __float128, __float128))
{
    return sqrtq (powq (x - p_x, 2) + powq (y - p_y, 2));
}

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
         double control (double, double, double, double, double, double))
{
    return cos (mult * x) * pow (control (t, x, y, p_x, p_y, mult), 2) -
           pow (y, 2);
}

static double
    test_x_d (double t, double x, double y, double p_x, double p_y, double mult,
              double control (double, double, double, double, double, double))
{
    return -y;
}

static double
    test_y_d (double t, double x, double y, double p_x, double p_y, double mult,
              double control (double, double, double, double, double, double))
{
    return x;
}

static double test_px_d (double t, double x, double y, double p_x, double p_y,
                         double mult,
                         double control (double, double, double, double, double,
                                         double))
{
    return -p_y;
}

static double test_py_d (double t, double x, double y, double p_x, double p_y,
                         double mult,
                         double control (double, double, double, double, double,
                                         double))
{
    return p_x;
}

static double
    test_B_d (double t, double x, double y, double p_x, double p_y, double mult,
              double control (double, double, double, double, double, double))
{
    return sqrt (pow (x - p_x, 2) + pow (y - p_y, 2));
}

#pragma GCC diagnostic pop
