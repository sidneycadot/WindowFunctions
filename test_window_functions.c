
/////////////////////////////
// test_window_functions.c //
/////////////////////////////

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "window_functions.h"

#define MAX_WINDOW_NAME_SIZE 100

static double * init_window(const char * window_name, unsigned n)
{
    double * w = NULL;

    if (strcmp(window_name, "rectwin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        rectwin(w, n);
    }
    else if (strcmp(window_name, "hann") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        hann(w, n, true);
    }
    else if (strcmp(window_name, "hann_periodic") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        hann(w, n, false);
    }
    else if (strcmp(window_name, "hann_symmetric") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        hann(w, n, true);
    }
    else if (strcmp(window_name, "hamming") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        hamming(w, n, true);
    }
    else if (strcmp(window_name, "hamming_periodic") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        hamming(w, n, false);
    }
    else if (strcmp(window_name, "hamming_symmetric") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        hamming(w, n, true);
    }
    else if (strcmp(window_name, "blackman") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        blackman(w, n, true);
    }
    else if (strcmp(window_name, "blackman_periodic") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        blackman(w, n, false);
    }
    else if (strcmp(window_name, "blackman_symmetric") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        blackman(w, n, true);
    }
    else if (strcmp(window_name, "blackmanharris") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        blackmanharris(w, n, true);
    }
    else if (strcmp(window_name, "blackmanharris_periodic") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        blackmanharris(w, n, false);
    }
    else if (strcmp(window_name, "blackmanharris_symmetric") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        blackmanharris(w, n, true);
    }
    else if (strcmp(window_name, "nuttallwin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        nuttallwin(w, n, true);
    }
    else if (strcmp(window_name, "nuttallwin_periodic") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        nuttallwin(w, n, false);
    }
    else if (strcmp(window_name, "nuttallwin_symmetric") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        nuttallwin(w, n, true);
    }
    else if (strcmp(window_name, "flattopwin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        flattopwin(w, n, true);
    }
    else if (strcmp(window_name, "flattopwin_periodic") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        flattopwin(w, n, false);
    }
    else if (strcmp(window_name, "flattopwin_symmetric") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        flattopwin(w, n, true);
    }
    else if (strcmp(window_name, "triang") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        triang(w, n);
    }
    else if (strcmp(window_name, "bartlett") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        bartlett(w, n);
    }
    else if (strcmp(window_name, "barthannwin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        barthannwin(w, n);
    }
    else if (strcmp(window_name, "bohmanwin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        bohmanwin(w, n);
    }
    else if (strcmp(window_name, "parzenwin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        parzenwin(w, n);
    }
    else if (strcmp(window_name, "gausswin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        gausswin(w, n, 2.5);
    }
    else if (strcmp(window_name, "gausswin_2p5") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        gausswin(w, n, 2.5);
    }
    else if (strcmp(window_name, "gausswin_3p2") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        gausswin(w, n, 3.2);
    }
    else if (strcmp(window_name, "tukeywin") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        tukeywin(w, n, 0.5);
    }
    else if (strcmp(window_name, "tukeywin_0p0") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        tukeywin(w, n, 0.0);
    }
    else if (strcmp(window_name, "tukeywin_0p2") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        tukeywin(w, n, 0.2);
    }
    else if (strcmp(window_name, "tukeywin_0p5") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        tukeywin(w, n, 0.5);
    }
    else if (strcmp(window_name, "tukeywin_0p8") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        tukeywin(w, n, 0.8);
    }
    else if (strcmp(window_name, "tukeywin_1p0") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        tukeywin(w, n, 1.0);
    }
    else if (strcmp(window_name, "kaiser") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        kaiser(w, n, 0.5);
    }
    else if (strcmp(window_name, "kaiser_0p5") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        kaiser(w, n, 0.5);
    }
    else if (strcmp(window_name, "kaiser_0p8") == 0)
    {
        w = (double *)malloc(sizeof(double) * n);
        kaiser(w, n, 0.8);
    }

    return w;
}

int main()
{
    const char * filename = "reference/matlab_windows.txt";

    FILE * f = fopen(filename, "r");

    char     current_window_name[MAX_WINDOW_NAME_SIZE] = "";
    unsigned current_window_n                          = 0;
    double * current_window                            = NULL;

    for (;;)
    {
        char window_name[MAX_WINDOW_NAME_SIZE];
        unsigned n;
        unsigned i;
        double reference_value;

        int result = fscanf(f, "%s%u%u%lf", window_name, &n, &i, &reference_value);
        if (result != 4)
        {
            break;
        }

        bool current_window_ok = (current_window != NULL && strcmp(current_window_name, window_name) == 0 && current_window_n == n);

        if (!current_window_ok)
        {
            double * calculated_window = init_window(window_name, n);
            if (calculated_window == NULL)
            {
                printf("WARNING: cannot init window \"%s\" with n == %u\n", window_name, n);
            }
            else
            {
                strcpy(current_window_name, window_name);

                free(current_window);

                current_window    = calculated_window;
                current_window_n  = n;
                current_window_ok = true;
            }
        }

        if (current_window_ok)
        {
            const double calculated_value = current_window[i - 1];

            const double error = calculated_value - reference_value;
            if (fabs(error) > 1e-10)
            {
                printf("ERROR: window \"%s\" n = %u i = %u : reference = %f, calculated = %f, error = %f\n", window_name, n, i, reference_value, calculated_value, error);
            }
        }
    }

    free(current_window);

    fclose(f);
}
