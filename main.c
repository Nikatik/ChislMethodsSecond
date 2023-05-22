#include "Functions.h"
#include "Includes.h"
#include "Newton.h"
#include "RungeKutta.h"

#include <time.h>

int    mains (int);
int    open (FILE**, char*, unsigned int*, unsigned int*);
void   initialization (long long unsigned, long long unsigned, unsigned int,
                       unsigned int, unsigned int, double***, double***, double**,
                       double**, double**);
void   reading (FILE*, double**, unsigned int, long long unsigned,
                long long unsigned);
double rd (FILE*);
void printing (double**, unsigned int, long long unsigned, long long unsigned);
void cleaning (long long unsigned, unsigned int, FILE*, double**, double**,
               double*, double*, double*);

int  mainq (int);
void initializationq (long long unsigned, long long unsigned, unsigned int,
                      unsigned int, unsigned int, __float128***, __float128***,
                      __float128**, __float128**, __float128**);
void readingq (FILE*, __float128**, unsigned int, long long unsigned,
               long long unsigned);
__float128 rdq (FILE*);
void       printingq (__float128**, unsigned int, long long unsigned,
                      long long unsigned);
void       cleaningq (long long unsigned, unsigned int, FILE*, __float128**,
                      __float128**, __float128*, __float128*, __float128*);

int main (int argc, char* argv[])
{
    if (argc > 1)
    {
        int n;
        sscanf (argv[1], "%d", &n);
        if (n < -5 && n >= -17) return mains (n);
        else if (n < -17) return mainq (n);
    }
    return mains (-12);
}

int mains (int toller)
{
    double tol = pow (10, toller);

    double x  = 1;
    double y  = 0;
    double px = 1;
    double py = 0;
    double b  = 0;

    FILE* inpf;

    double*      T;
    unsigned int T_c = 1;

    double*      mult;
    unsigned int mult_c = 3;

    double* alpha;

    unsigned int s;
    unsigned int p;

    double** k;
    double** cab;

    clock_t timer;

    long long unsigned i = 0;
    long long unsigned j = 0;

    /////////////////////////////////////////////////////////////////////////////////////////

    switch (open (&inpf, "koef (8).txt", &s, &p))
    {
        case -1: return -1;
        case -2: return -2;
        default: printf ("File opened.\n");
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    initialization (i, j, s, T_c, mult_c, &cab, &k, &T, &mult, &alpha);
    printf ("All initialized.\n");

    /////////////////////////////////////////////////////////////////////////////////////////

    reading (inpf, cab, s, i, j);
    printf ("RK matrix readed.\n\n");

    // printing (cab, s, i, j);

    /////////////////////////////////////////////////////////////////////////////////////////

    // Testing ARK
    astep (pi * pow (10, 3), &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol,
           test_x_d, test_y_d, test_px_d, test_py_d, test_B_d, 1, 1, 0);
    printf (
        "The Runge-Kutta test:\n%.2e    %.2e  |  %.2e    %.2e  |  %.7e  |  10^3*Pi\n\n",
        x - 1, y, px - 1, py, b);

    /////////////////////////////////////////////////////////////////////////////////////////

    // The shoting method
    timer = clock();

    x  = 11;
    y  = 0;
    px = 0;
    py = 0;

    for (unsigned int m = 0; m < T_c; m++)
        for (unsigned int l = 0; l < mult_c; l++)
        {
            alpha[4] = mult[l];
            alpha[5] = T[m];

            if (shooting_method (500, T[m], x, py, alpha, p, s, k, cab, tol,
                                 x_d, y_d, px_d, py_d, B_d, mult[l]) == -1)
                break;
        }

    timer -= clock();
    printf ("%.5f seconds\n", ((double) -timer) / CLOCKS_PER_SEC);

    /////////////////////////////////////////////////////////////////////////////////////////

    cleaning (i, s, inpf, cab, k, T, mult, alpha);
    printf ("All cleaned.\n");

    /////////////////////////////////////////////////////////////////////////////////////////

    return 0;
}

// File checker
int open (FILE** inpf, char* path, unsigned int* s, unsigned int* p)
{

    *inpf = fopen (path, "r");
    if (*inpf == NULL)
    {
        printf ("File doen`t exist\n");
        return -1;
    }
    if (!fscanf (*inpf, "%ud", p) || !fscanf (*inpf, "%ud", s))
    {
        printf ("File isn`t correct\n");
        fclose (*inpf);
        return -2;
    }
    return 0;
}

// Constructor
void initialization (long long unsigned i, long long unsigned j, unsigned int s,
                     unsigned int T_c, unsigned int mult_c, double*** cab,
                     double*** k, double** T, double** mult, double** alpha)
{
    *cab = (double**) malloc ((s + 3) * sizeof (double*));
    for (i = 0; i < s + 3; i++)
    {
        (*cab)[i] = (double*) malloc (s * sizeof (double));
        for (j = 0; j < s; j++) (*cab)[i][j] = 0;
    }


    *k = (double**) malloc (5 * sizeof (double*));
    for (i = 0; i < 5; i++)
    {
        (*k)[i] = (double*) malloc (s * sizeof (double));
        for (j = 0; j < s; j++) (*k)[i][j] = 0;
    }

    // values for T from book
    *T = (double*) malloc (T_c * sizeof (double));
    for (i = 0; i < T_c; i++)
    {
        switch (i)
        {
            case 0: (*T)[i] = 4; break;
            case 1: (*T)[i] = 2; break;
            case 2: (*T)[i] = 3; break;
            case 3: (*T)[i] = 3.5; break;
            case 4: (*T)[i] = 10; break;
            case 5: (*T)[i] = 0.1; break;
            case 7: (*T)[i] = 0.1; break;
            default: (*T)[i] = 1;
        }
    }

    // values for Alpha from book
    *mult = (double*) malloc (mult_c * sizeof (double));
    for (i = 0; i < mult_c; i++)
    {
        switch (i)
        {
            case 0: (*mult)[i] = 0.0; break;
            default:
                (*mult)[i] = (i <= 4) ? pow (10, (double) i - 4.)
                                      : 5. * ((double) i - 4.);
        }
    }

    *alpha      = (double*) malloc (6 * sizeof (double));
    (*alpha)[0] = (*mult)[0];
    (*alpha)[1] = (*T)[0];
    (*alpha)[2] = -9.;        // start value for shoting method
    (*alpha)[3] = 3.;         // start value for shoting method
}

// CAB matrix reading
void reading (FILE* inpf, double** cab, unsigned int s, long long unsigned i,
              long long unsigned j)
{
    for (i = 0; i < s; i++) cab[0][i] = rd (inpf);

    for (i = 2; i < s + 3; i++)
        for (j = 0; j + 1 < i && j < s && !feof (inpf); j++)
            cab[i][j] = rd (inpf);
}

// Reading floating point number from a/b format
double rd (FILE* inpf)
{
    int  chisl = 0;
    int  znam  = 1;
    char tmp;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"

    fscanf (inpf, "%d", &chisl);
    if (fscanf (inpf, "%c", &tmp) && tmp == '/') fscanf (inpf, "%d", &znam);

#pragma GCC diagnostic pop

    return ((double) chisl) / znam;
}

// CAB matrix printing
void printing (double** cab, unsigned int s, long long unsigned i,
               long long unsigned j)
{
    printf ("\n");
    for (i = 0; i < s + 3; i++)
    {
        for (j = 0; j < s; j++) printf ("%6.3f ", cab[i][j]);
        printf ("\n");
    }
    printf ("\n");
}

// Deconstructor
void cleaning (long long unsigned i, unsigned int s, FILE* inpf, double** cab,
               double** k, double* T, double* mult, double* alpha)
{
    fclose (inpf);

    for (i = 0; i < s + 3; i++) free (cab[i]);
    free (cab);

    for (i = 0; i < 5; i++) free (k[i]);
    free (k);

    free (T);
    free (mult);
    free (alpha);
}

int mainq (int toller)
{
    __float128 tol = (__float128) pow (10, toller);

    __float128 x  = 1;
    __float128 y  = 0;
    __float128 px = 1;
    __float128 py = 0;
    __float128 b  = 0;

    FILE* inpf;

    __float128*  T;
    unsigned int T_c = 1;

    __float128*  mult;
    unsigned int mult_c = 3;

    __float128* alpha;

    unsigned int s;
    unsigned int p;

    __float128** k;
    __float128** cab;

    clock_t timer;

    long long unsigned i = 0;
    long long unsigned j = 0;

    /////////////////////////////////////////////////////////////////////////////////////////

    switch (open (&inpf, "koef (8).txt", &s, &p))
    {
        case -1: return -1;
        case -2: return -2;
        default: printf ("File opened.\n");
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    initializationq (i, j, s, T_c, mult_c, &cab, &k, &T, &mult, &alpha);
    printf ("All initialized.\n");

    /////////////////////////////////////////////////////////////////////////////////////////

    readingq (inpf, cab, s, i, j);
    printf ("RK matrix readed.\n\n");

    // printing (cab, s, i, j);

    /////////////////////////////////////////////////////////////////////////////////////////

    // Testing ARK
    astepq (M_PIq * (__float128) pow (10, 3), &x, &y, &px, &py, &b, &i, &j, p,
            s, k, cab, tol, test_x_dq, test_y_dq, test_px_dq, test_py_dq,
            test_B_dq, 1, 1, 0);
    printf (
        "The Runge-Kutta test:\n%.2Le    %.2Le  |  %.2Le    %.2Le  |  %.7Le  |  10^3*Pi\n\n",
        (long double) x - 1, (long double) y, (long double) px - 1,
        (long double) py, (long double) b);

    /////////////////////////////////////////////////////////////////////////////////////////

    // The shoting method
    timer = clock();

    x  = 11;
    y  = 0;
    px = 0;
    py = 0;

    for (unsigned int m = 0; m < T_c; m++)
        for (unsigned int l = 0; l < mult_c; l++)
        {
            alpha[4] = mult[l];
            alpha[5] = T[m];

            if (shooting_methodq (500, T[m], x, py, alpha, p, s, k, cab, tol,
                                  x_dq, y_dq, px_dq, py_dq, B_dq,
                                  mult[l]) == -1)
                break;
        }

    timer -= clock();
    printf ("%.5Lf seconds\n", ((long double) -timer) / CLOCKS_PER_SEC);

    /////////////////////////////////////////////////////////////////////////////////////////

    cleaningq (i, s, inpf, cab, k, T, mult, alpha);
    printf ("All cleaned.\n");

    /////////////////////////////////////////////////////////////////////////////////////////

    return 0;
}

// Constructor
void initializationq (long long unsigned i, long long unsigned j,
                      unsigned int s, unsigned int T_c, unsigned int mult_c,
                      __float128*** cab, __float128*** k, __float128** T,
                      __float128** mult, __float128** alpha)
{
    *cab = (__float128**) malloc ((s + 3) * sizeof (__float128*));
    for (i = 0; i < s + 3; i++)
    {
        (*cab)[i] = (__float128*) malloc (s * sizeof (__float128));
        for (j = 0; j < s; j++) (*cab)[i][j] = 0;
    }


    *k = (__float128**) malloc (5 * sizeof (__float128*));
    for (i = 0; i < 5; i++)
    {
        (*k)[i] = (__float128*) malloc (s * sizeof (__float128));
        for (j = 0; j < s; j++) (*k)[i][j] = 0;
    }

    // values for T from book
    *T = (__float128*) malloc (T_c * sizeof (__float128));
    for (i = 0; i < T_c; i++)
    {
        switch (i)
        {
            case 0: (*T)[i] = 4; break;
            case 1: (*T)[i] = 2; break;
            case 2: (*T)[i] = 3; break;
            case 3: (*T)[i] = 3.5Q; break;
            case 4: (*T)[i] = 10; break;
            case 5: (*T)[i] = 0.1Q; break;
            case 7: (*T)[i] = 0.1Q; break;
            default: (*T)[i] = 1;
        }
    }

    // values for Alpha from book
    *mult = (__float128*) malloc (mult_c * sizeof (__float128));
    for (i = 0; i < mult_c; i++)
    {
        switch (i)
        {
            case 0: (*mult)[i] = 0; break;
            default:
                (*mult)[i] = (i <= 4) ? powq (10.Q, i - 4.Q)
                                      : 5.Q * (i - 4.Q);
        }
        //printf("%Le\n",(long double) (*mult)[i]);
    }

    *alpha      = (__float128*) malloc (6 * sizeof (__float128));
    (*alpha)[0] = (*mult)[0];
    (*alpha)[1] = (*T)[0];
    (*alpha)[2] = -9;        // start value for shoting method
    (*alpha)[3] = 3;         // start value for shoting method
}

// CAB matrix reading
void readingq (FILE* inpf, __float128** cab, unsigned int s,
               long long unsigned i, long long unsigned j)
{
    for (i = 0; i < s; i++) cab[0][i] = rdq (inpf);

    for (i = 2; i < s + 3; i++)
        for (j = 0; j + 1 < i && j < s && !feof (inpf); j++)
            cab[i][j] = rdq (inpf);
}

// Reading floating point number from a/b format
__float128 rdq (FILE* inpf)
{
    int  chisl = 0;
    int  znam  = 1;
    char tmp;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"

    fscanf (inpf, "%d", &chisl);
    if (fscanf (inpf, "%c", &tmp) && tmp == '/') fscanf (inpf, "%d", &znam);

#pragma GCC diagnostic pop

    return ((__float128) chisl) / znam;
}

// CAB matrix printing
void printingq (__float128** cab, unsigned int s, long long unsigned i,
                long long unsigned j)
{
    printf ("\n");
    for (i = 0; i < s + 3; i++)
    {
        for (j = 0; j < s; j++) printf ("%6.3Lf ", (long double) cab[i][j]);
        printf ("\n");
    }
    printf ("\n");
}

// Deconstructor
void cleaningq (long long unsigned i, unsigned int s, FILE* inpf,
                __float128** cab, __float128** k, __float128* T,
                __float128* mult, __float128* alpha)
{
    fclose (inpf);

    for (i = 0; i < s + 3; i++) free (cab[i]);
    free (cab);

    for (i = 0; i < 5; i++) free (k[i]);
    free (k);

    free (T);
    free (mult);
    free (alpha);
}
