
////////////////
// test_fft.c //
////////////////

#include <stdio.h>
#include <complex.h>
#include <math.h>

#include "window_functions.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

static void forward_transform(complex double * z, unsigned size)
{
    complex double zcopy[size];

    for (unsigned i = 0; i < size; ++i)
    {
        zcopy[i] = z[i];
    }

    for (unsigned i = 0; i < size; ++i)
    {
        z[i] = 0;
        for (unsigned j = 0; j < size; ++j)
        {
            z[i] += cexp(-2 * M_PI * I * i * j / size) * zcopy[j];
        }
    }
}

static void inverse_transform(complex double * z, unsigned size)
{
    complex double zcopy[size];

    for (unsigned i = 0; i < size; ++i)
    {
        zcopy[i] = z[i];
    }

    for (unsigned i = 0; i < size; ++i)
    {
        z[i] = 0;
        for (unsigned j = 0; j < size; ++j)
        {
            z[i] += cexp(+2 * M_PI * I * i * j / size) / size * zcopy[j];
        }
    }
}

static void init_test_array(complex double * z, unsigned size)
{
    for (unsigned i = 0; i < size; ++i)
    {
        double re_numerator   = (i % 11) - 5;
        double re_denominator = (i % 13) + 1;
        double im_numerator   = (i % 17) - 8;
        double im_denominator = (i % 19) + 1;

        z[i] = (re_numerator / re_denominator) + (im_numerator / im_denominator) * I;
    }
}

int main()
{
    for (unsigned size = 1; size <= 16; ++size)
    {
        complex double z1[size];
        complex double z2[size];

        // Test forward FFT

        init_test_array(z1, size);
        init_test_array(z2, size);

        forward_transform(z1, size);
        fft((double *)z2, size, false);

        double forward_fft_max_err = 0;

        for (unsigned i = 0; i < size; ++i)
        {
            const double err = cabs(z2[i] - z1[i]);

            forward_fft_max_err = fmax(forward_fft_max_err, err);
        }

        // Test inverse FFT

        init_test_array(z1, size);
        init_test_array(z2, size);

        inverse_transform(z1, size);
        fft((double *)z2, size, true);

        double inverse_fft_max_err = 0;

        for (unsigned i = 0; i < size; ++i)
        {
            const double err = cabs(z2[i] - z1[i]);

            inverse_fft_max_err = fmax(inverse_fft_max_err, err);
        }

        // Present results

        printf("size: %6u    forward_fft_max_err: %12.6g    inverse_fft_max_err: %12.6g\n", size,
               forward_fft_max_err, inverse_fft_max_err);
    }

    return 0;
}
