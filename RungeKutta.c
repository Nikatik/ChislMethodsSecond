#include "RungeKutta.h"

#include "Control.h"

// Runge-Kutta in general form (one step)
void RK (double t, double* x0, double* y0, double* px0, double* py0, double* b0,
         double h, unsigned int s, double** k, double** cab,
         double f (double, double, double, double, double, double,
                   double (*) (double, double, double, double, double, double)),
         double g (double, double, double, double, double, double,
                   double (*) (double, double, double, double, double, double)),
         double u (double, double, double, double, double, double,
                   double (*) (double, double, double, double, double, double)),
         double v (double, double, double, double, double, double,
                   double (*) (double, double, double, double, double, double)),
         double l (double, double, double, double, double, double,
                   double (*) (double, double, double, double, double, double)),
         double mult, double c (double, double, double, double, double, double),
         bool check)
{
    double temp_x;
    double temp_y;
    double temp_px;
    double temp_py;
    double temp_b;


    for (unsigned int i = 0; i < s; i++)
    {
        temp_x  = 0;
        temp_y  = 0;
        temp_px = 0;
        temp_py = 0;
        temp_b  = 0;
        for (unsigned int j = 0; j < i; j++)
        {
            temp_x += cab[i + 1][j] * k[0][j];
            temp_y += cab[i + 1][j] * k[1][j];
            temp_px += cab[i + 1][j] * k[2][j];
            temp_py += cab[i + 1][j] * k[3][j];
            temp_b += cab[i + 1][j] * k[4][j];
        }
        // printf("%f\n",f(1,0));
        k[0][i] = f (t + h * cab[0][i], *x0 + h * temp_x, *y0 + h * temp_y,
                     *px0 + h * temp_px, *py0 + h * temp_py, mult, c);

        k[1][i] = g (t + h * cab[0][i], *x0 + h * temp_x, *y0 + h * temp_y,
                     *px0 + h * temp_px, *py0 + h * temp_py, mult, c);

        k[2][i] = u (t + h * cab[0][i], *x0 + h * temp_x, *y0 + h * temp_y,
                     *px0 + h * temp_px, *py0 + h * temp_py, mult, c);

        k[3][i] = v (t + h * cab[0][i], *x0 + h * temp_x, *y0 + h * temp_y,
                     *px0 + h * temp_px, *py0 + h * temp_py, mult, c);

        k[4][i] = l (t + h * cab[0][i], *x0 + h * temp_x, *y0 + h * temp_y,
                     *px0 + h * temp_px, *py0 + h * temp_py, mult, c);
    }

    temp_x  = 0;
    temp_y  = 0;
    temp_px = 0;
    temp_py = 0;
    temp_b  = 0;


    if (check)
    {
        for (unsigned int i = 0; i < s; i++)
        {
            temp_x += cab[s + 2][i] * k[0][i];
            temp_y += cab[s + 2][i] * k[1][i];
            temp_px += cab[s + 2][i] * k[2][i];
            temp_py += cab[s + 2][i] * k[3][i];
            temp_b += cab[s + 2][i] * k[4][i];
        }
    }
    else
    {
        for (unsigned int i = 0; i < s; i++)
        {
            temp_x += cab[s + 1][i] * k[0][i];
            temp_y += cab[s + 1][i] * k[1][i];
            temp_px += cab[s + 1][i] * k[2][i];
            temp_py += cab[s + 1][i] * k[3][i];
            temp_b += cab[s + 1][i] * k[4][i];
        }
    }

    *x0 += h * temp_x;
    *y0 += h * temp_y;
    *px0 += h * temp_px;
    *py0 += h * temp_py;
    *b0 += h * temp_b;
}

// Finding crossing point
void diagonal (
    double dist, double x, double y, double px, double py, double b,
    double* temp_x, double* temp_y, double* temp_px, double* temp_py,
    double* temp_b, double* temp, unsigned int s, double** k, double** cab,
    double tol,
    double f (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double g (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double u (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double v (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double l (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double mult, double c (double, double, double, double, double, double),
    double func (double, double, double, double, double, double), bool print)
{
    double lx = x, ly = y, lpx = px, lpy = py, lb = b, ldist = dist;
    double _x = x, _y = y, _px = px, _py = py, _b = b, _dist = dist, _h = *temp;

    for (; fabs (func (_dist, _x, _y, _px, _py, mult)) >
           tol * pow (10, 5);)        // finding point of crossing
    {
        if (func (_dist, _x, _y, _px, _py, mult) * func (dist + *temp, *temp_x,
                                                         *temp_y, *temp_px,
                                                         *temp_py, mult) <
            0)
        {
            lx    = _x;
            ly    = _y;
            lpx   = _px;
            lpy   = _py;
            lb    = _b;
            ldist = _dist;
            _h    = fabs (func (_dist, _x, _y, _px, _py, mult)) /
                 fabs (func (dist + *temp, *temp_x, *temp_y, *temp_px, *temp_py,
                             mult) -
                       func (_dist, _x, _y, _px, _py, mult)) *
                 (dist + *temp - _dist);        // step from _x to startx

            RK (_dist, &_x, &_y, &_px, &_py, &_b, _h, s, k, cab, f, g, u, v, l,
                mult, *c, false);
            _dist += _h;
        }
        else
        {
            *temp_x  = _x;
            *temp_y  = _y;
            *temp_px = _px;
            *temp_py = _py;
            *temp_b  = _b;
            *temp    = _dist - dist;

            _h = fabs (func (ldist, lx, ly, lpx, lpy, mult)) /
                 fabs (func (ldist, lx, ly, lpx, lpy, mult) -
                       func (_dist, _x, _y, _px, _py, mult)) *
                 (_dist - ldist);        // step from lx to startx

            _x    = lx;
            _y    = ly;
            _px   = lpx;
            _py   = lpy;
            _b    = lb;
            _dist = ldist;

            RK (_dist, &_x, &_y, &_px, &_py, &_b, _h, s, k, cab, f, g, u, v, l,
                mult, *c, false);
            _dist += _h;
        }
    }

    *temp_x  = x;
    *temp_y  = y;
    *temp_px = px;
    *temp_py = py;
    *temp_b  = b;
    *temp    = _dist - dist;

    RK (dist, temp_x, temp_y, temp_px, temp_py, temp_b, *temp, s, k, cab, f, g,
        u, v, l, mult, *c, false);        // counting startx point
    if (print)
        printf (
            "Turning control point:\tt = %8.4f,\tX = %8.4f,\tY = %8.4f,\tPx = %8.4f,\tPy = %8.4f,\tB = %8.4f,\tDiscrepancy: %13e\n",
            dist + *temp, *temp_x, *temp_y, *temp_px, *temp_py, *temp_b,
            func (dist + *temp, *temp_x, *temp_y, *temp_px, *temp_py, mult));
}

// Adaptive Runge-Kutta
double astep (
    double T, double* x, double* y, double* px, double* py, double* b,

    long long unsigned* i, long long unsigned* j,

    unsigned int p, unsigned int s,

    double** k, double** cab,

    double tol,

    double f (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double g (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double u (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double v (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double l (double, double, double, double, double, double,
              double (*) (double, double, double, double, double, double)),
    double mult, bool test, bool print)
{
    double temp_x, temp_y, temp_px, temp_py, temp_b, x_, y_, px_, py_, b_, dist,
        h, temp, fac, err, norm;
    double (*c) (double, double, double, double, double, double);
    x_ = temp_x = *x;
    y_ = temp_y = *y;

    px_ = temp_px = *px;
    py_ = temp_py = *py;
    b_ = temp_b = *b;


    dist = 0;
    h    = 0.01;
    fac  = 1.7;
    err  = 0;
    norm = 0;
    c    = (test)                  ? 0
           : (cos (mult * *x) > 0) ? control_n
           : (*py > 0)             ? control_p
                                   : control_m;
    // printf("%.2e    %.6e  |  %.6e    %.2e\n", *x, *y, *px, *py);

    for (*i = *j = 0; T - dist > EPS;)
    {

        if (dist + h > T) h = T - dist;

        RK (dist, &temp_x, &temp_y, &temp_px, &temp_py, &temp_b, h, s, k, cab,
            f, g, u, v, l, mult, *c, false);
        RK (dist, &x_, &y_, &px_, &py_, &b_, h, s, k, cab, f, g, u, v, l, mult,
            *c, true);


        // norm = fmax( fmax( fabs (temp_x - x_), fabs (temp_y - y_)), fmax(fabs
        // (temp_px - px_), fabs (temp_py - py_)));
        norm = sqrt (pow (temp_x - x_, 2) + pow (temp_y - y_, 2) +
                     pow (temp_px - px_, 2) + pow (temp_py - py_, 2));
        temp = h;
        h *= fmin (fac, fmax (0.7, pow (0.98 * tol / norm, 1. / (p + 1))));
        if (h < pow (10, -18))
        {
            // printf ("\nSomething goes wrong...\n");
            return -1000000;
        }
        if (norm > tol)
        {
            x_ = temp_x = *x;
            y_ = temp_y = *y;

            px_ = temp_px = *px;
            py_ = temp_py = *py;
            b_ = temp_b = *b;
            fac         = 1;
            *j += 1;
            continue;
        }

        if (fabs (cos (mult * *x)) > tol * pow (10, 6) &&
            cos (mult * *x) * cos (mult * temp_x) < 0)
        {
            diagonal (dist, *x, *y, *px, *py, *b, &temp_x, &temp_y, &temp_px,
                      &temp_py, &temp_b, &temp, s, k, cab, tol, f, g, u, v, l,
                      mult, *c, turn_x, print);

            x_  = *x;
            y_  = *y;
            px_ = *px;
            py_ = *py;
            b_  = *b;

            RK (dist, &x_, &y_, &px_, &py_, &b_, temp, s, k, cab, f, g, u, v, l,
                mult, *c, true);
            norm = sqrt (pow (temp_x - x_, 2) + pow (temp_y - y_, 2) +
                         pow (temp_px - px_, 2) + pow (temp_py - py_, 2));
            c    = (temp_py > 0) ? control_p : control_m;
        }
        if (fabs (cos (mult * *x)) < tol * pow (10, 6) &&
            fabs (cos (mult * temp_x)) > tol * pow (10, 6))
            c = (cos (mult * temp_x) > 0) ? control_n
                : (temp_py > 0)           ? control_p
                                          : control_m;

        if (fabs (*py) > tol * pow (10, 6) && *py * temp_py < 0)
        {
            diagonal (dist, *x, *y, *px, *py, *b, &temp_x, &temp_y, &temp_px,
                      &temp_py, &temp_b, &temp, s, k, cab, tol, f, g, u, v, l,
                      mult, *c, turn_py, print);

            x_  = *x;
            y_  = *y;
            px_ = *px;
            py_ = *py;
            b_  = *b;

            RK (dist, &x_, &y_, &px_, &py_, &b_, temp, s, k, cab, f, g, u, v, l,
                mult, *c, true);
            norm = sqrt (pow (temp_x - x_, 2) + pow (temp_y - y_, 2) +
                         pow (temp_px - px_, 2) + pow (temp_py - py_, 2));
            c    = (cos (mult * temp_x) > 0) ? control_n : control_p;
        }
        if (fabs (*py) < tol * pow (10, 6) &&
            fabs (temp_py) > tol * pow (10, 6))
            c = (cos (mult * temp_x) > 0) ? control_n
                : (temp_py > 0)           ? control_p
                                          : control_m;

        err += norm;
        dist += temp;
        *x  = temp_x;
        *y  = temp_y;
        *px = temp_px;
        *py = temp_py;
        *b  = temp_b;
        fac = 1.7;
        *i += 1;
    }
    // printf("%.2e    %.2e  |  %.2e    %.2e\n", *x, *y, *px, *py);
    return err;
}
