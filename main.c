#include <time.h>

#include "Includes.h"
#include "Functions.h"
#include "Newton.h"
#include "RungeKutta.h"

    
double rd (FILE* inpf)        // reading floating point number from a/b format
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

int main()
{
    double tol           = pow (10, -17);
    long long unsigned i = 0;
    long long unsigned j = 0;

    double x=1;
	double y=0;
	double px=1;
	double py=0;
    double b=0;
	
    double*     T; // T из условия
    unsigned int T_c = 3;
    double*  mult; // alpha из условия
    unsigned int mult_c = 1;
    double* alpha; // под невязки

    unsigned int s;
    unsigned int p;
    double** k;
    double** cab;

    clock_t timer;

    /////////////////////////////////////////////////////////////////////////////////////////

    // file check


    FILE* inpf = fopen ("koef (8).txt", "r");
    if (inpf == NULL)
    {
        printf ("File doen`t exist\n");
        return -1;
    }
    if (!fscanf (inpf, "%ud", &p) || !fscanf (inpf, "%ud", &s))
    {
        printf ("File isn`t correct\n");
        fclose (inpf);
        return -2;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    // CAB matrix initialization

    cab = (double**)malloc ((s + 3) * sizeof (double*));
    for (i = 0; i < s + 3; i++)
    {
        cab[i] = (double*)malloc (s * sizeof (double));
        for (j = 0; j < s; j++)
        {
            cab[i][j] = 0;
        }
    }
    
    T = (double*)malloc (T_c * sizeof (double));
    for(i=0; i < T_c; i++){
        switch (i) {
            case 0:
                T[i] = 1;
                break;
            case 1:
                T[i] = 2;
                break;
            case 2:
                T[i] = 3;
                break;
            case 3:
                T[i] = 3.5;
                break;
            case 4:
                T[i] = 10;
                break;
            case 5:
                T[i] = 0.1;
                break;
            case 7:
                T[i] = 0.1;
                break;
            default:
                T[i] = 1;
        }
    }
    
    mult = (double*)malloc (mult_c * sizeof (double));
    for(i=0; i < mult_c; i++){
        switch (i) {
            case 0:
                mult[i] = 0.0;
                break;
            default:
                mult[i] = (i<=4)?pow(10,(double)i - 4.):5.*((double)i - 4.);
        }
    }

	alpha = (double*)malloc (6 * sizeof (double));
    alpha[0] = mult[0];
    alpha[1] = T[0];
    alpha[2] = 5.;
    alpha[3] = -12.;

    k    = (double**)malloc (5 * sizeof (double*));
    k[0] = (double*)malloc (s * sizeof (double));
    k[1] = (double*)malloc (s * sizeof (double));
	k[2] = (double*)malloc (s * sizeof (double));
    k[3] = (double*)malloc (s * sizeof (double));
    k[4] = (double*)malloc (s * sizeof (double));
	
    for (i = 0; i < s; i++)
    {
        k[0][i] = 0;
        k[1][i] = 0;
		k[2][i] = 0;
        k[3][i] = 0;
        k[4][i] = 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    // CAB matrix reading

    for (i = 0; i < s; i++)
    {
        cab[0][i] = rd (inpf);
    }

    for (i = 2; i < s + 3; i++)
        for (j = 0; j + 1 < i && j < s && !feof (inpf); j++)
            cab[i][j] = rd (inpf);

    // CAB matrix printing
    /*
            printf ("\n");
            for (unsigned int i = 0; i < s + 3; i++)
            {
                for (unsigned int j = 0; j < s; j++)
                    printf ("%6.3f ", cab[i][j]);
                printf ("\n");
            }
            printf ("\n");
        */
    ///////////////////////////////////////////////////////////////////////////////////////// 
    
    // the shooting method
	astep(pi*pow(10,3), &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol, test_x_d, test_y_d, test_px_d, test_py_d, test_B_d, 1, 1, 0);
    printf("The Runge-Kutta test:\n%.2e    %.2e  |  %.2e    %.2e  |  %.7e  |  10^3*Pi\n\n", x - 1, y, px - 1, py, b);
    
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

    /////////////////////////////////////////////////////////////////////////////////////////

    // cleaning :)

    fclose (inpf);
    for (i = 0; i < s + 3; i++)
    {
        free (cab[i]);
    }
    free (cab);
	free (T);
	free (mult);
	free (alpha);
    free (k[0]);
    free (k[1]);
	free (k[2]);
    free (k[3]);
    free (k[4]);
    free (k);

    timer -= clock();
    printf ("%.5f seconds\n", ((double)-timer)/CLOCKS_PER_SEC);
    return 0;
}
