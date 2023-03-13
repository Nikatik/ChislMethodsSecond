#include <math.h>
#include <stdio.h>

_Pragma ("GCC diagnostic push")
_Pragma ("GCC diagnostic ignored \"-Wunused-parameter\"")
static double contol_n (double t, double x, double y, double p_x, double p_y, double mult)
{
    if(p_y > 0)
        return (p_y/cos(mult*x)<24)?p_y/cos(mult*x):24;
    else
        return (p_y/cos(mult*x)>-24)?p_y/cos(mult*x):-24;
}
/*
static double contol_p (double t, double x, double y, double p_x, double p_y, double mult)
{
    return 24;
}

static double contol_m (double t, double x, double y, double p_x, double p_y, double mult)
{
    return -24;
}
*/
_Pragma ("GCC diagnostic pop")