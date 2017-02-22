
/////////////////////////////
// time_window_functions.c //
/////////////////////////////

#include <stdio.h>
#include <time.h>

#include "window_functions.h"

int main()
{
    unsigned size = 4096;

    double w[size];
    clock_t t1 = clock();
    chebwin(w, size, 100.0);
    clock_t t2 = clock();

    clock_t duration = (t2 - t1);

    double duration_ms = 1000.0 * duration / CLOCKS_PER_SEC;

    printf("duration: %10.3f ms\n", duration_ms);

    return 0;
}
