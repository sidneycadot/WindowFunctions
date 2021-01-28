# About the WindowFunctions

This repository provides a C implementation for generating signal processing windows.
It reproduces the windows provided in Matlab and Octave up to numerical precision.

As a bonus, a complex-valued, in-place FFT implementation is included that works for
any size *n* in O(*n* log *n*) time.

## Design decisions

The implementation is written in C99. It can also be used from C++, if you first compile
the C file using a C99-compliant compiler.

It is assumed that math functions are available, including functions that
implement complex math. All functions used are defined in the C99 standard.

Calculations of Chebyshev windows depend on the ability to perform Fourier
transforms, and calculation of Kaiser windows depend on the availability of
the modified Bessel function I0(x). Implementations of these functions are
included, which means that there are no external dependencies.

The implementation does not use dynamic stack allocation (malloc), enhancing
suitability for application in embedded applications. However, if you intend
to use the 'chebwin()' function, see its specific notes below, since it may
temporarily use a considerable amount of stack space.

All floating-point calculations are performed using double precision.

## Performance

The implementation was optimized for correctness and simplicity of implementation
rather than performance. The exception is that we do provide an FFT implementation
in order to support the chebwin() function. (An implementation using a direct
Fourier transform would have been simpler, but unbearably slow for all but the
smallest sizes of windows.)

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
(Gauss, Tukey, Taylor, Kaiser, Chebyshev) have tunable parameters. The corresponding
functions take extra arguments in the function call.

## List of implemented windows

### Cosine Windows

Cosine windows consist of a linear combination of cosine-shaped base-functions. The domain
of the window is x = [0 .. 2*pi], and the base functions are { cos(0 * x), cos(1 * x),
cos(2 * x), cos (3 * x), ... }, etcetera. Note that the first function cos (0 * x) simplifies to a
constant 1.

The number of cosine terms used in a cosine window is referred to as its *order*.

The cosine window functions have an *sflag* parameter that specifies whether the generated
window should be 'symmetric' (sflag == true), or 'periodic' (sflag == false).
The symmetric choice returns a perfectly symmetric window.
The periodic choice works as if a window with (size + 1) elements is calculated,
after which the last element is dropped.

#### Generic cosine window 

The library provides a function to synthesize a generic cosine window by specifying a list
of coefficients:

    void cosine_window(double * w, unsigned n, const double * coeff, unsigned ncoeff, bool sflag);

#### Specific cosine windows

The first-order 'rectwin' function gives a rectangular or uniform window, which can be seen as a degenerate
case of a cosine window: w = cos(0 * x). Note that this function lacks the 'sflag' parameter
since it would make no difference.

    void rectwin(double * w, unsigned n);

The Hann window is a second-order cosine window that is shaped as a single period of the cosine function.

    void hann(double * w, unsigned n, bool sflag);

The Hamming windows is a second-order cosine window that is shaped as a cosine like the Hann window, except it
tapers off to the edge value 0.08 rather than 0:

    void hamming(double * w, unsigned n, bool sflag);

The Blackman window is a third-order cosine window:

    void blackman(double * w, unsigned n, bool sflag);

The Blackman-Harris window is a fourth-order cosine window:

    void blackmanharris(double * w, unsigned n, bool sflag);

The Nuttall window is a fourth-order window that is very similar to the Blackman-Harris window.
The Nutall window definitions are different between Matlab and Octave. We define both:

    void nuttallwin(double * w, unsigned n, bool sflag);
    void nuttallwin_octave(double * w, unsigned n, bool sflag);

The signal processing literature defines several 'flattop' windows. Matlab and Octave define similar
but not idential fifth-order cosine windows with this name; we define both.

These windows are characterized by a very flat top, and by the fact that they contain negative values:

    void flattopwin(double * w, unsigned n, bool sflag);
    void flattopwin_octave(double * w, unsigned n, bool sflag);

### OTHER WINDOWS, NOT PARAMETERIZED

void triang            (double * w, unsigned n);
void bartlett          (double * w, unsigned n);
void barthannwin       (double * w, unsigned n);
void bohmanwin         (double * w, unsigned n);
void parzenwin         (double * w, unsigned n);

### OTHER WINDOWS, PARAMETERIZED

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
from a complex double pointer to a double pointer on invocation.
