#include "Newton.h"

#include "RungeKutta.h"

// Discrepancy counting
double error (
    double T, double x, double y, double px, double py, double b, double* err,
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
    double mult, bool print)
{
    long long unsigned i = 0;
    long long unsigned j = 0;

    if (astep (T, &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol, f, g, u, v,
               l, mult, 0, print) +
            1000000 <
        EPS)
        return -1000000;

    err[0] = x;
    err[1] = y - 1;

    return b;
}

// The shooting method
int shooting_method (
    unsigned int max_iterations, double T, double x, double py, double* alpha,
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
    double mult)
{
    double alpha_1;        // x
    double alpha_2;        // py
    double betta;          // dempher

    unsigned int iteration;

    double** A = (double**) malloc (2 * sizeof (double*));
    A[0]       = (double*) malloc (2 * sizeof (double));
    A[0][0]    = 0;
    A[0][1]    = 0;
    A[1]       = (double*) malloc (2 * sizeof (double));
    A[1][0]    = 0;
    A[1][1]    = 0;

    double* B = (double*) malloc (6 * sizeof (double));
    for (int i = 0; i < 6; i++) B[i] = 0;

    double* newton = (double*) malloc (2 * sizeof (double));
    newton[0]      = 0;
    newton[1]      = 0;

    double* err = (double*) malloc (2 * sizeof (double));
    err[0]      = 0;
    err[1]      = 0;

    double delta  = pow (2, -11);
    double step_m = pow (10, -1) * ((alpha[0] - alpha[4] > 0) ? -1 : 1);
    double step_T = pow (10, -1) * ((alpha[1] - alpha[5] > 0) ? -1 : 1);
    double integral;

    for (; 1;)
    {
        integral = 0;
        alpha_1  = alpha[2];        // y
        alpha_2 =
            alpha[3];        // px
                             //  printf("mult_0=%13.7e  |  mult_k=%13.7e  |
                             //  mult=%13.7e  |  T_0=%13.7e  |  T_k=%13.7e  |
                             //  T=%13.7e  |  step_m=%13.7e  |
                             //  step_T=%13.7e\n",alpha[0], alpha[4], mult,
                             //  alpha[1], alpha[5], T, step_m, step_T);

        for (iteration = 0;; iteration++)
        {
            // printf("T=%14.7e  |  Y(0)=%14.7e  |  Px(0)=%14.7e  |
            // Mult=%14.7e\n",alpha[5], alpha_1, alpha_2, alpha[4]);
            if (iteration == max_iterations)
            {
                if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS) step_m *= 0.1;
                if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS) step_T *= 0.1;
                break;
            }
            if (iteration % 5 == 0)
            {
                // dy

                if (error (alpha[5], x, alpha_1 + delta, alpha_2, py, integral,
                           err, p, s, k, cab, tol, f, g, u, v, l, alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[0] = err[0];        // y
                B[3] = err[1];        // px

                if (error (alpha[5], x, alpha_1 - delta, alpha_2, py, integral,
                           err, p, s, k, cab, tol, f, g, u, v, l, alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[1] = err[0];        // y
                B[4] = err[1];        // px

                if (error (alpha[5], x, alpha_1 + 2 * delta, alpha_2, py,
                           integral, err, p, s, k, cab, tol, f, g, u, v, l,
                           alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[2] = err[0];        // y
                B[5] = err[1];        // px

                if (error (alpha[5], x, alpha_1 - 2 * delta, alpha_2, py,
                           integral, err, p, s, k, cab, tol, f, g, u, v, l,
                           alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                A[0][0] = (2. * B[0] - 2. * B[1] - B[2] / 4. + err[0] / 4.) /
                          (3 * delta);
                A[1][0] = (2. * B[3] - 2. * B[4] - B[5] / 4. + err[1] / 4.) /
                          (3 * delta);

                // dpx

                if (error (alpha[5], x, alpha_1, alpha_2 + delta, py, integral,
                           err, p, s, k, cab, tol, f, g, u, v, l, alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[0] = err[0];        // y
                B[3] = err[1];        // px

                if (error (alpha[5], x, alpha_1, alpha_2 - delta, py, integral,
                           err, p, s, k, cab, tol, f, g, u, v, l, alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[1] = err[0];        // y
                B[4] = err[1];        // px

                if (error (alpha[5], x, alpha_1, alpha_2 + 2. * delta, py,
                           integral, err, p, s, k, cab, tol, f, g, u, v, l,
                           alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[2] = err[0];        // y
                B[5] = err[1];        // px

                if (error (alpha[5], x, alpha_1, alpha_2 - 2. * delta, py,
                           integral, err, p, s, k, cab, tol, f, g, u, v, l,
                           alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                A[0][1] = (2. * B[0] - 2. * B[1] - B[2] / 4. + err[0] / 4.) /
                          (3 * delta);
                A[1][1] = (2. * B[3] - 2. * B[4] - B[5] / 4. + err[1] / 4.) /
                          (3 * delta);
            }

            integral = 0;
            integral = error (alpha[5], x, alpha_1, alpha_2, py, integral, err,
                              p, s, k, cab, tol, f, g, u, v, l, alpha[4], 0);
            if (integral + 1000000 < EPS)
            {
                if (fabs (alpha[4] - (alpha[0] + step_m)) < EPS) step_m *= 0.1;
                if (fabs (alpha[5] - (alpha[1] + step_T)) < EPS) step_T *= 0.1;
                break;
            }

            B[0] = err[0];        // y
            B[1] = err[1];        // px

            // printf("%10.2e    %10.2e  |  %10.2e\n%10.2e    %10.2e  |
            // %10.2e\n", A[0][0], A[0][1], B[0], A[1][0], A[1][1], B[1]);

            B[5] = sqrt (pow (B[0], 2) + pow (B[1], 2));
            if (B[5] < tol * pow (10, 3))
            {
                // printf and break
                // printf("Shooting iteration: %u,\nAlpha: %.4f,\nDistance:
                // %.4f,\nY(0) = %.14f,\nPx(0) = %.14f,\nB =
                // %.14f,\nDiscrepancy: %.3e\n\n", iteration, alpha[4],
                // alpha[5], alpha_1, alpha_2, integral, sqrt(pow(B[0], 2) +
                // pow(B[1], 2)));

                alpha[0] = alpha[4];
                alpha[1] = alpha[5];
                alpha[2] = alpha_1;
                alpha[3] = alpha_2;

                break;
            }

            // printf("Shooting iteration: %ud,\nAlpha: %.3f,\nY(0) =
            // %.3e,\nPx(0) = %.3e,\nDiscrepancy: %.3e\n\n", iteration,
            // mult,alpha_1,alpha_2,sqrt(pow(B[0],2)+pow(B[1],2)));

            newton[0] = (A[0][0] * B[0] - A[1][0] * B[1]) /
                        (A[1][0] * A[0][1] - A[0][0] * A[1][1]);
            newton[1] = (A[0][0] * B[1] - A[1][0] * B[0]) /
                        (A[1][0] * A[0][1] - A[0][0] * A[1][1]);

            // demph

            for (betta = 1;;)
            {
                if (error (alpha[5], x, alpha_1 - betta * newton[0],
                           alpha_2 + betta * newton[1], py, integral, err, p, s,
                           k, cab, tol, f, g, u, v, l, alpha[4], 0) +
                        1000000 <
                    EPS)
                {
                    betta *= 1.01;
                    break;
                }
                if (sqrt (pow (err[0], 2) + pow (err[1], 2)) - B[5] < EPS &&
                    betta > pow (2, -15))
                {
                    B[5] = sqrt (pow (err[0], 2) + pow (err[1], 2));
                    betta /= 1.01;
                }
                else
                {
                    betta *= 1.01;
                    break;
                }
            }


            alpha_1 -= betta * newton[0];
            alpha_2 += betta * newton[1];
        }

        if (fabs (alpha[0] - alpha[4]) < EPS &&
            fabs (alpha[1] - alpha[5]) < EPS)
        {
            step_m *= 1.1;
            step_T *= 1.1;
            if (mult - alpha[4] < step_m) step_m = mult - alpha[4];

            if (T - alpha[5] < step_T) step_T = T - alpha[5];
        }

        if (fabs (mult - alpha[0]) > EPS) alpha[4] = alpha[0] + step_m;
        if (fabs (T - alpha[1]) > EPS) alpha[5] = alpha[1] + step_T;


        if (fabs (alpha[0] - mult) < EPS && fabs (alpha[1] - T) < EPS)
        {
            printf ("\n-----------------------\n");
            integral = 0;
            integral = error (alpha[5], x, alpha_1, alpha_2, py, integral, err,
                              p, s, k, cab, tol, f, g, u, v, l, alpha[4], 1);
            printf (
                "Shooting iteration: %u,\nAlpha: %.4f,\nDistance: %.4f,\nY(0) = %.14f,\nPx(0) = %.14f,\nB = %.14f,\nDiscrepancy: %.3e\n\n",
                iteration + 1, alpha[4], alpha[5], alpha_1, alpha_2, integral,
                sqrt (pow (B[0], 2) + pow (B[1], 2)));
            break;
        }

        if ((fabs (alpha[5] - alpha[1]) < pow (10, -4) &&
             fabs (alpha[5] - alpha[1]) > EPS) ||
            (fabs (alpha[4] - alpha[0]) < pow (10, -4) &&
             fabs (alpha[4] - alpha[0]) > EPS))
        {
            printf ("\n-----------------------------------\n");
            printf ("Newton's method does not converge!\n");
            integral = 0;
            integral =
                error (alpha[1], x, alpha[2], alpha[3], py, integral, err, p, s,
                       k, cab, tol, f, g, u, v, l, alpha[0], 1);
            printf (
                "Shooting iteration: %u,\nAlpha: %16.4f,\nRequested alpha: %.4f,\nDistance: %16.4f,\nRequested distance: %.4f,\nY(0) = %.14f,\nPx(0) = %.14f,\nB = %.14f,\nDiscrepancy: %.3e\nSteps: %.4e   %.4e\n\n",
                iteration + 1, alpha[0], mult, alpha[1], T, alpha[2], alpha[3],
                integral, sqrt (pow (err[0], 2) + pow (err[1], 2)), step_m,
                step_T);

            free (A[0]);
            free (A[1]);
            free (A);
            free (newton);
            free (B);
            free (err);

            return -1;
        }
        iteration = 0;
    }

    free (A[0]);
    free (A[1]);
    free (A);
    free (newton);
    free (B);
    free (err);

    return 0;
}