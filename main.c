#include <time.h>

#include "Includes.h"
#include "Functions.h"
#include "Newton.h"
#include "RungeKutta.h"

int open (FILE**, char*, unsigned int*, unsigned int*);
void initialization (long long unsigned, long long unsigned, unsigned int, unsigned int, unsigned int, double***, double***, double**, double**, double**);
void reading (FILE*, double**, unsigned int, long long unsigned, long long unsigned);
double rd (FILE*);
void printing (double**, unsigned int, long long unsigned, long long unsigned);
void cleaning (long long unsigned, unsigned int, FILE*, double**, double**, double*, double*, double*);

int main()
{
    double tol = pow (10, -17);

    double x  = 1;
	double y  = 0;
	double px = 1;
	double py = 0;
    double b  = 0;
    
    FILE* inpf;
	
    double*      T;
    unsigned int T_c = 3;
    
    double*      mult;
    unsigned int mult_c = 1;
    
    double* alpha;

    unsigned int s;
    unsigned int p;
    
    double** k;
    double** cab;
    
    clock_t timer;
    
    long long unsigned i = 0;
    long long unsigned j = 0;

    /////////////////////////////////////////////////////////////////////////////////////////

    switch (open(&inpf, "koef (8).txt", &s, &p)) 
    {
        case -1:
            return -1;
        case -2:
            return -2;
        default:;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    initialization (i, j, s, T_c, mult_c, &cab, &k, &T, &mult, &alpha);
    
    /////////////////////////////////////////////////////////////////////////////////////////

    reading(inpf, cab, s, i, j);

    // printing (cab, s, i, j);
    
    ///////////////////////////////////////////////////////////////////////////////////////// 
    
    // Testing ARK
	astep(pi*pow(10,3), &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol, test_x_d, test_y_d, test_px_d, test_py_d, test_B_d, 1, 1, 0);
    printf("The Runge-Kutta test:\n%.2e    %.2e  |  %.2e    %.2e  |  %.7e  |  10^3*Pi\n\n", x - 1, y, px - 1, py, b);
    
    ///////////////////////////////////////////////////////////////////////////////////////// 
    
    // The shoting method
    timer = clock();
    
	x=0;
	y=0;
	px=0;
	py=0;
    
    for(unsigned int m = 0; m < T_c; m++)
        for(unsigned int l = 0; l < mult_c; l++){
            alpha[4] = mult[l];
            alpha[5] = T[m];
        
            if(shooting_method(500, T[m], x, py, alpha, p, s, k, cab, tol, x_d, y_d, px_d, py_d, B_d, mult[l]) == -1)
                break;
        }
    
    timer -= clock();
    printf ("%.5f seconds\n", ((double)-timer)/CLOCKS_PER_SEC);

    /////////////////////////////////////////////////////////////////////////////////////////

    cleaning(i, s, inpf, cab, k, T, mult, alpha);
    
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
void initialization (long long unsigned i, long long unsigned j, unsigned int s, unsigned int T_c, unsigned int mult_c, double*** cab, double*** k, double** T, double** mult, double** alpha)
{
    *cab = (double**)malloc ((s + 3) * sizeof (double*));
    for (i = 0; i < s + 3; i++)
    {
        (*cab)[i] = (double*)malloc (s * sizeof (double));
        for (j = 0; j < s; j++)
        {
            (*cab)[i][j] = 0;
        }
    }
    
    
    *k = (double**)malloc (5 * sizeof (double*));
    for (i = 0; i < 5; i++)
    {
        (*k)[i] = (double*)malloc (s * sizeof (double));
        for (j = 0; j < s; j++)
        {
            (*k)[i][j] = 0;
        }
    }
    
    //values for T from book
    *T = (double*)malloc (T_c * sizeof (double));
    for(i=0; i < T_c; i++){
        switch (i) {
            case 0:
                (*T)[i] = 1;
                break;
            case 1:
                (*T)[i] = 2;
                break;
            case 2:
                (*T)[i] = 3;
                break;
            case 3:
                (*T)[i] = 3.5;
                break;
            case 4:
                (*T)[i] = 10;
                break;
            case 5:
                (*T)[i] = 0.1;
                break;
            case 7:
                (*T)[i] = 0.1;
                break;
            default:
                (*T)[i] = 1;
        }
    }
    
    //values for Alpha from book
    *mult = (double*)malloc (mult_c * sizeof (double));
    for(i=0; i < mult_c; i++){
        switch (i) {
            case 0:
                (*mult)[i] = 0.0;
                break;
            default:
                (*mult)[i] = (i<=4)?pow(10,(double)i - 4.):5.*((double)i - 4.);
        }
    }

	*alpha = (double*)malloc (6 * sizeof (double));
    (*alpha)[0] = (*mult)[0];
    (*alpha)[1] = (*T)[0];
    (*alpha)[2] = -0.5; //start value for shoting method
    (*alpha)[3] = -2.;  //start value for shoting method
}

// CAB matrix reading
void reading (FILE* inpf, double** cab, unsigned int s, long long unsigned i, long long unsigned j)
{
    for (i = 0; i < s; i++)
    {
        cab[0][i] = rd (inpf);
    }

    for (i = 2; i < s + 3; i++)
        for (j = 0; j + 1 < i && j < s && !feof (inpf); j++)
            cab[i][j] = rd (inpf);
}

// Reading floating point number from a/b format
double rd (FILE* inpf)        
{
    int chisl = 0;
    int znam  = 1;
    char tmp;
    _Pragma ("GCC diagnostic push");
    _Pragma ("GCC diagnostic ignored \"-Wunused-result\"");
    fscanf (inpf, "%d", &chisl);
    if (fscanf (inpf, "%c", &tmp) && tmp == '/') fscanf (inpf, "%d", &znam);
    _Pragma ("GCC diagnostic pop");
    return ((double)chisl) / znam;
}

// CAB matrix printing
void printing (double** cab, unsigned int s, long long unsigned i, long long unsigned j)
{
    printf ("\n");
    for (i = 0; i < s + 3; i++)
    {
        for (j = 0; j < s; j++)
            printf ("%6.3f ", cab[i][j]);
        printf ("\n");
    }
    printf ("\n");
}

// Deconstructor
void cleaning (long long unsigned i, unsigned int s, FILE* inpf, double** cab, double** k, double* T, double* mult, double* alpha)
{
    fclose (inpf);
    
    for (i = 0; i < s + 3; i++)
        free (cab[i]);
    free (cab);
    
    for (i = 0; i < 5; i++)
        free (k[i]);
    free (k);
    
	free (T);
	free (mult);
	free (alpha);
}
