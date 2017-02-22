# About the WindowFunctions

This repository provides C functions for the generation of signal processing windows.
It reproduces the windows as found in Matlab and Octave up to numerical precision.

As a bonus, a complex, in-place FFT implementation is provided that works for any
size *n* in O(*n* log *n*) time.

## Design decisions

The implementation is written in C99. It can also be used from C++.

It is assumed that math functions are available, including functions that
implement complex math. All functions used are defined in the C99 standard.

Calculations of Chebyshev windows depend on the ability to perform Fourier
transforms, and calculation of Kaiser windows depend on the availability of
the modified Bessel function I0(x). Implementation of these functions is
included, which means that there are no external dependencies.

The implementation does not use dynamic stack allocation (malloc), enhancing
suitability for application in embedded applications. However, if you intend to
use the 'chebwin()' function, see its specific notes below, since it may
temporarily use a considerable amount of stack space.

All floating-point calculations are performed using double precision.

## Performance

The implementation was optimized for correctness and simplicity of implementation
rather than performance. The exception is that we do provide an FFT implementation
in order to support the chebwin() function. (An implementation using a direct
Fourier transform would have been simpler, but unbearably slow for all but the smallest
sizes of windows.)

On a Cortex-M4 class ARM processor lacking double precision floating point support,
running at 100 MHz, generation of a 4096-point Chebyshev window takes about 3 seconds.
On a multi-GHz Intel processor, this same operation takes 3 milliseconds.

## Using the window functions

The caller of any of the window functions is responsible for reserving memory for storage
of the window values.

All window functions take a pointer-to-double as their first argument, which defines where
the calculated window values will be stored. The second argument 'size' defines the size
of the window, i.e., the number of elements.

Most window types are fully defined by their size parameter. However, several window types
(Gauss, Tukey, Taylor, Kaiser, Chebyshev) have tunable parameters. These are provided as
extra arguments in the function call.

## List of implemented windows

// COSINE WINDOWS

// Generic cosine window

void cosine_window     (double * w, unsigned n, const double * coeff, unsigned ncoeff, bool sflag);

// Specific cosine windows

void rectwin           (double * w, unsigned n            ); // order == 1
void hann              (double * w, unsigned n, bool sflag); // order == 2
void hamming           (double * w, unsigned n, bool sflag); // order == 3
void blackman          (double * w, unsigned n, bool sflag); // order == 3
void blackmanharris    (double * w, unsigned n, bool sflag); // order == 4
void nuttallwin        (double * w, unsigned n, bool sflag); // order == 4
void nuttallwin_octave (double * w, unsigned n, bool sflag); // order == 4
void flattopwin        (double * w, unsigned n, bool sflag); // order == 5
void flattopwin_octave (double * w, unsigned n, bool sflag); // order == 5

// OTHER WINDOWS, NOT PARAMETERIZED

void triang            (double * w, unsigned n);
void bartlett          (double * w, unsigned n);
void barthannwin       (double * w, unsigned n);
void bohmanwin         (double * w, unsigned n);
void parzenwin         (double * w, unsigned n);

// OTHER WINDOWS, PARAMETERIZED

void gausswin          (double * w, unsigned n, double alpha);
void tukeywin          (double * w, unsigned n, double r);
void taylorwin         (double * w, unsigned n, unsigned nbar, double sll);
void kaiser            (double * w, unsigned n, double beta);
void chebwin           (double * w, unsigned n, double r);

### Chebyshev window

The only window function that uses any type of allocation is the Chebyshev window function
'chebwin'. It allocates memory on the stack using C99's Variable-length array (VLA) feature. See
the notes below for more information.

## Tests

TBW

## FFT

Calculation of Chebyshev windows depends on the ability to perform Fourier transforms. Since
this functionality had to be implemented anyway, is non-trivial, and potentially useful, it is
made available as a function.

The FFT function takes three parameters: a pointer to an array containing (2 * size) doubles,
a parameter 'size', and a boolean flag 'inv' determining whether a forward (inv == false) or
inverse (inv == true) Discrete Fourier Transform is to be performed.

Note that the array pointer is a pointer-to-double rather than a pointer-to-complex-double.
This enables the same prototype to be used from C++, which defines complex numbers via
templates rather than as a built-in type in the language. Use casting to convert to type
from a complex double pointer to a double pointer on ionvocation.
