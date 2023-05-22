#include <math.h>
#include <quadmath.h>
#include <stdio.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

static __float128 control_nq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
                         __float128 mult)
{
    if (p_y > 0) return (p_y / cosq (mult * x) < 24) ? p_y / cosq (mult * x) : 24;
    else return (p_y / cosq (mult * x) > -24) ? p_y / cosq (mult * x) : -24;
}

static __float128 control_pq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
                         __float128 mult)
{
    return 24;
}

static __float128 control_mq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
                         __float128 mult)
{
    return -24;
}

static __float128 turn_xq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
                      __float128 mult)
{
    return cosq (mult * x);
}

static __float128 turn_pyq (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
                       __float128 mult)
{
    return p_y;
}

/*
static __float128 hamelt (__float128 t, __float128 x, __float128 y, __float128 p_x, __float128 p_y,
__float128 mult, __float128 c)
{
    return -cos(mult * x) * pow(c, 2.) / 2. + p_y * c;
}
*/

static double control_n (double t, double x, double y, double p_x, double p_y,
                         double mult)
{
    if (p_y > 0) return (p_y / cos (mult * x) < 24) ? p_y / cos (mult * x) : 24;
    else return (p_y / cos (mult * x) > -24) ? p_y / cos (mult * x) : -24;
}

static double control_p (double t, double x, double y, double p_x, double p_y,
                         double mult)
{
    return 24;
}

static double control_m (double t, double x, double y, double p_x, double p_y,
                         double mult)
{
    return -24;
}

static double turn_x (double t, double x, double y, double p_x, double p_y,
                      double mult)
{
    return cos (mult * x);
}

static double turn_py (double t, double x, double y, double p_x, double p_y,
                       double mult)
{
    return p_y;
}

/*
static double hamelt (double t, double x, double y, double p_x, double p_y,
double mult, double c)
{
    return -cos(mult * x) * pow(c, 2.) / 2. + p_y * c;
}
*/
#pragma GCC diagnostic pop
