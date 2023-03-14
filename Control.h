#include <math.h>
#include <stdio.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

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

/*
static double hamelt (double t, double x, double y, double p_x, double p_y, double mult, double c)
{
    return -cos(mult * x) * pow(c, 2.) / 2. + p_y * c;
}
*/
#pragma GCC diagnostic pop
