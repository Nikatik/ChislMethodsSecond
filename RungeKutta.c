#include "RungeKutta.h"

#include "Control.h"

// Runge-Kutta in general form (one step)
void RKq (__float128 t, __float128* x0, __float128* y0, __float128* px0, __float128* py0, __float128* b0,
         __float128 h, unsigned int s, __float128** k, __float128** cab,
         __float128 f (__float128, __float128, __float128, __float128, __float128, __float128,
                   __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
         __float128 g (__float128, __float128, __float128, __float128, __float128, __float128,
                   __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
         __float128 u (__float128, __float128, __float128, __float128, __float128, __float128,
                   __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
         __float128 v (__float128, __float128, __float128, __float128, __float128, __float128,
                   __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
         __float128 l (__float128, __float128, __float128, __float128, __float128, __float128,
                   __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
         __float128 mult, __float128 c (__float128, __float128, __float128, __float128, __float128, __float128),
         bool check)
{
    __float128 temp_x;
    __float128 temp_y;
    __float128 temp_px;
    __float128 temp_py;
    __float128 temp_b;


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
void diagonalq (
    __float128 dist, __float128 x, __float128 y, __float128 px, __float128 py, __float128 b,
    __float128* temp_x, __float128* temp_y, __float128* temp_px, __float128* temp_py,
    __float128* temp_b, __float128* temp, unsigned int s, __float128** k, __float128** cab,
    __float128 tol,
    __float128 f (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 g (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 u (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 v (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 l (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 mult, __float128 c (__float128, __float128, __float128, __float128, __float128, __float128),
    __float128 func (__float128, __float128, __float128, __float128, __float128, __float128), bool print)
{
    __float128 lx = x, ly = y, lpx = px, lpy = py, lb = b, ldist = dist;
    __float128 _x = x, _y = y, _px = px, _py = py, _b = b, _dist = dist, _h = *temp;

    for (; fabsq (func (_dist, _x, _y, _px, _py, mult)) >
           tol * powq (tol, -1.Q / 3.Q);)        // finding point of crossing
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
            _h    = fabsq (func (_dist, _x, _y, _px, _py, mult)) /
                 fabsq (func (dist + *temp, *temp_x, *temp_y, *temp_px, *temp_py,
                             mult) -
                       func (_dist, _x, _y, _px, _py, mult)) *
                 (dist + *temp - _dist);        // step from _x to startx

            RKq (_dist, &_x, &_y, &_px, &_py, &_b, _h, s, k, cab, f, g, u, v, l,
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

            _h = fabsq (func (ldist, lx, ly, lpx, lpy, mult)) /
                 fabsq (func (ldist, lx, ly, lpx, lpy, mult) -
                       func (_dist, _x, _y, _px, _py, mult)) *
                 (_dist - ldist);        // step from lx to startx

            _x    = lx;
            _y    = ly;
            _px   = lpx;
            _py   = lpy;
            _b    = lb;
            _dist = ldist;

            RKq (_dist, &_x, &_y, &_px, &_py, &_b, _h, s, k, cab, f, g, u, v, l,
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

    RKq (dist, temp_x, temp_y, temp_px, temp_py, temp_b, *temp, s, k, cab, f, g,
        u, v, l, mult, *c, false);        // counting startx point
    if (print)
        printf (
            "Turning control point:\tt = %8.4Lf,\tX = %8.4Lf,\tY = %8.4Lf,\tPx = %8.4Lf,\tPy = %8.4Lf,\tB = %8.4Lf,\tDiscrepancy: %13Le\n",
            (long double)(dist + *temp), (long double)*temp_x, (long double)*temp_y, (long double)*temp_px, (long double)*temp_py, (long double)*temp_b,
            (long double)func (dist + *temp, *temp_x, *temp_y, *temp_px, *temp_py, mult));
}

// Adaptive Runge-Kutta
__float128 astepq (
    __float128 T, __float128* x, __float128* y, __float128* px, __float128* py, __float128* b,

    long long unsigned* i, long long unsigned* j,

    unsigned int p, unsigned int s,

    __float128** k, __float128** cab,

    __float128 tol,

    __float128 f (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 g (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 u (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 v (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 l (__float128, __float128, __float128, __float128, __float128, __float128,
              __float128 (*) (__float128, __float128, __float128, __float128, __float128, __float128)),
    __float128 mult, bool test, bool print)
{
    __float128 temp_x, temp_y, temp_px, temp_py, temp_b, x_, y_, px_, py_, b_, dist,
        h, temp, fac, err, norm;
    __float128 (*c) (__float128, __float128, __float128, __float128, __float128, __float128);
    x_ = temp_x = *x;
    y_ = temp_y = *y;

    px_ = temp_px = *px;
    py_ = temp_py = *py;
    b_ = temp_b = *b;


    dist = 0;
    h    = 0.01Q;
    fac  = 1.7Q;
    err  = 0;
    norm = 0;
    c    = (test)                  ? 0
           : (cosq (mult * *x) > 0) ? control_nq
           : (*py > 0)             ? control_pq
                                   : control_mq;
    // printf("%.2e    %.6e  |  %.6e    %.2e\n", *x, *y, *px, *py);

    for (*i = *j = 0; T - dist > EPSq;)
    {

        if (dist + h > T) h = T - dist;

        RKq (dist, &temp_x, &temp_y, &temp_px, &temp_py, &temp_b, h, s, k, cab,
            f, g, u, v, l, mult, *c, false);
        RKq (dist, &x_, &y_, &px_, &py_, &b_, h, s, k, cab, f, g, u, v, l, mult,
            *c, true);


        // norm = fmax( fmax( fabs (temp_x - x_), fabs (temp_y - y_)), fmax(fabs
        // (temp_px - px_), fabs (temp_py - py_)));
        norm = sqrtq (powq (temp_x - x_, 2) + powq (temp_y - y_, 2) +
                     powq (temp_px - px_, 2) + powq (temp_py - py_, 2));
        temp = h;
        h *= fminq (fac, fmaxq (0.7Q, powq (0.98Q * tol / norm, 1.Q / (p + 1))));
        if (h < (__float128) pow (10, -17))
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

        if (fabsq (cosq (mult * *x)) > tol * powq (tol, -12.Q / 30.Q) &&
            cosq (mult * *x) * cosq (mult * temp_x) < 0)
        {
            diagonalq (dist, *x, *y, *px, *py, *b, &temp_x, &temp_y, &temp_px,
                      &temp_py, &temp_b, &temp, s, k, cab, tol, f, g, u, v, l,
                      mult, *c, turn_xq, print);

            x_  = *x;
            y_  = *y;
            px_ = *px;
            py_ = *py;
            b_  = *b;

            RKq (dist, &x_, &y_, &px_, &py_, &b_, temp, s, k, cab, f, g, u, v, l,
                mult, *c, true);
            norm = sqrtq (powq (temp_x - x_, 2) + powq (temp_y - y_, 2) +
                         powq (temp_px - px_, 2) + powq (temp_py - py_, 2));
            c    = (temp_py > 0) ? control_pq : control_mq;
        }
        if (fabsq (cosq (mult * *x)) < tol * powq (tol, -12.Q / 30.Q) &&
            fabsq (cosq (mult * temp_x)) > tol * powq (tol, -12.Q / 30.Q))
            c = (cosq (mult * temp_x) > 0) ? control_nq
                : (temp_py > 0)           ? control_pq
                                          : control_mq;

        if (fabsq (*py) > tol * powq (tol, -12.Q / 30.Q) && *py * temp_py < 0)
        {
            diagonalq (dist, *x, *y, *px, *py, *b, &temp_x, &temp_y, &temp_px,
                      &temp_py, &temp_b, &temp, s, k, cab, tol, f, g, u, v, l,
                      mult, *c, turn_pyq, print);

            x_  = *x;
            y_  = *y;
            px_ = *px;
            py_ = *py;
            b_  = *b;

            RKq (dist, &x_, &y_, &px_, &py_, &b_, temp, s, k, cab, f, g, u, v, l,
                mult, *c, true);
            norm = sqrtq (powq (temp_x - x_, 2) + powq (temp_y - y_, 2) +
                         powq (temp_px - px_, 2) + powq (temp_py - py_, 2));
            c    = (cosq (mult * temp_x) > 0) ? control_nq : control_pq;
        }
        if (fabsq (*py) < tol * powq (tol, -12.Q / 30.Q) &&
            fabsq (temp_py) > tol * powq (tol, -12.Q / 30.Q))
            c = (cosq (mult * temp_x) > 0) ? control_nq
                : (temp_py > 0)           ? control_pq
                                          : control_mq;

        err += norm;
        dist += temp;
        *x  = temp_x;
        *y  = temp_y;
        *px = temp_px;
        *py = temp_py;
        *b  = temp_b;
        fac = 1.7Q;
        *i += 1;
    }
    // printf("%.2e    %.2e  |  %.2e    %.2e\n", *x, *y, *px, *py);
    return err;
}

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
           tol * pow (tol, -1. / 3.);)        // finding point of crossing
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
        if (h < pow (10, -17))
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

        if (fabs (cos (mult * *x)) > tol * pow (tol, -12. / 30.) &&
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
        if (fabs (cos (mult * *x)) < tol * pow (tol, -12. / 30.) &&
            fabs (cos (mult * temp_x)) > tol * pow (tol, -12. / 30.))
            c = (cos (mult * temp_x) > 0) ? control_n
                : (temp_py > 0)           ? control_p
                                          : control_m;

        if (fabs (*py) > tol * pow (tol, -12. / 30.) && *py * temp_py < 0)
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
        if (fabs (*py) < tol * pow (tol, -12. / 30.) &&
            fabs (temp_py) > tol * pow (tol, -12. / 30.))
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
