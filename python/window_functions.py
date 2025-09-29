"""Window functions implemented in Python with numpy.

These functions aim to reproduce the behavior of the Matlab functions with the same name.

The implementations below use only basic NumPy functions; the aim is to have reference implementations
for all window functions from which we can easily derive a C implementation.

There are two window types for which this is somewhat tricky:
- The Kaiser window depends on the ability to calculate the modified Bessel function of the first kind I0(x).
- The Chebyshev window depends on the ability to perform a Discrete Fourier Transform (DFT).
"""

import numpy as np

from bessel_i0 import bessel_i0  # Needed for the Kaiser window.

########################
#                      #
#    COSINE WINDOWS    #
#                      #
########################

def cosine_window(n: int, coefficients: list[float], symmetry_flag: bool=True) -> np.ndarray:
    """Generalized cosine window.

    Many window functions described in signal processing literature
    can be written as linear combinations of cosines over the window length.

    Let 'x' be values going from 0 for the first element, and 2*pi for the last element.

    The window can then be written as:

    w = c0 * cos(0 * x) + c1 * cos(1 * x) + c2 * cos(2 * x) + c3 * cos(3 * x) + ...

    (Note that the first term simplifies to just the constant value c0.)

    Examples of cosine windows implemented in Matlab:

                                    c0          c1           c2           c3            c4
        -------------------------------------------------------------------------------------------
        rectangular window          1.0
        hann window                 0.5         -0.5
        hamming window              0.54        -0.46
        blackman window             0.42        -0.5         0.08
        blackman-harris window      0.35875     -0.48829     0.14128      -0.01168
        nuttall window              0.3635819   -0.4891775   0.1365995    -0.0106411
        flattop window              0.21557895  -0.41663158  0.277263158  -0.083578947  0.006947368

        The "flattop" coefficients given above follow Matlab's "flattopwin" implementation.

        Note 1 : Octave defines the flattopwin coefficients differently, see implementation below.

                 The coefficient values used correspond to:

                 [0.21550795224343777, -0.4159303478298349, 0.2780052583940347, -0.08361708547045386, 0.006939356062238697]

        Note 2 : Octave defines the nuttallwin coefficients differently, see implementation below:

                 The coefficient values used are:

                 [0.355768, -0.487396, 0.144232, -0.012604]
    """

    assert isinstance(n, int) and n >= 1

    # Special case for (n==1), otherwise we'd divide by zero.
    if n == 1:
        return np.ones(1)

    if symmetry_flag:
        # The normal, symmetric form of the window.
        x = np.arange(n) * (2 * np.pi / (n - 1))
    else:
        # The periodic form. Calculate the window as if it's one element longer, and discard the last element.
        x = np.arange(n) * (2 * np.pi / (n    ))

    w = np.zeros(n)
    for (i, coefficient) in enumerate(coefficients):
        w += coefficient * np.cos(i * x)

    return w


def rectwin(n: int) -> np.ndarray:
    """Rectangular (all-ones) window.

    Note: this window can also be seen as the simplest of the cosine windows.
    It could be implemented as such:

        cosine_window(n, [1.0], True)

    Of course, the implementation given below is faster.
    """

    assert isinstance(n, int) and n >= 1

    return np.ones(n) # All values are 1. The dtype is 'float'.

def hann(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Hann window.

    Extrema are 0.
    Center value is 1 for odd length,
        0.5 - 0.5 * cos(pi * n / (n - 1)) for even length.
    """

    assert isinstance(n, int) and n >= 1

    coefficients = [0.5, -0.5]

    return cosine_window(n, coefficients, symmetry_flag)


def hanning(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Hanning window (Matlab version).

    In Matlab, the Hanning window is a truncated version of the Hann window.
    """

    assert isinstance(n, int) and n >= 1

    if n == 1:
        return np.ones(1)

    if symmetry_flag:
        return hann(n + 2)[1:-1].copy()
    else:
        return hann(n + 1)[0:-1].copy()


def hanning_octave(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Hanning window (Octave version).

    In Octave, the hanning() function is identical to the hann() function.
    """

    assert isinstance(n, int) and n >= 1

    coefficients = [0.5, -0.5]

    return cosine_window(n, coefficients, symmetry_flag)


def general_hamming_scipy(n: int, alpha: float, symmetry_flag: bool=True) -> np.ndarray:
    """Generalized Hamming window."""

    assert isinstance(n, int) and n >= 1

    coefficients = [alpha, alpha - 1.0]

    return cosine_window(n, coefficients, symmetry_flag)


def hamming(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Hamming window.

    Note that the Hamming window is raised; its extreme values are 0.08.

    The center value is 1 for odd length;
    The center values are 0.54 - 0.46 * cos(pi * n / (n - 1)) for even length.

    Note: Mathematica defines the Hamming window with slightly different coefficients:
      [25/46, 21/46].

    """

    assert isinstance(n, int) and n >= 1

    coefficients = [0.54, -0.46]

    return cosine_window(n, coefficients, symmetry_flag)


def blackman(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Blackman window."""

    assert isinstance(n, int) and n >= 1

    coefficients = [0.42, -0.5, 0.08]

    return cosine_window(n, coefficients, symmetry_flag)


def blackmanharris(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Blackman-Harris window.

    Note: this window is very similar to the Nuttall window.
    """

    assert isinstance(n, int) and n >= 1

    coefficients = [0.35875, -0.48829, 0.14128, -0.01168]

    return cosine_window(n, coefficients, symmetry_flag)


def nuttallwin(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Nuttall window (Matlab-compatible version).

    Note: this window is very similar to the Blackman-Harris window.
    Note: this window is available as the 'blackmannuttall' window in Octave.
    """

    assert isinstance(n, int) and n >= 1

    coefficients = [0.3635819, -0.4891775, 0.1365995, -0.0106411]

    return cosine_window(n, coefficients, symmetry_flag)


def nuttallwin_octave(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Nuttall window (Octave-compatible version)."""

    assert isinstance(n, int) and n >= 1

    coefficients = [0.355768, -0.487396, 0.144232, -0.012604]

    return cosine_window(n, coefficients, symmetry_flag)


def flattopwin(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Flattop window (Matlab-compatible version).

    Note: this window contains negative entries!
    """

    assert isinstance(n, int) and n >= 1

    coefficients = [0.21557895, -0.41663158, 0.277263158, -0.083578947, 0.006947368]

    return cosine_window(n, coefficients, symmetry_flag)


def flattopwin_octave(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Flattop window (Octave-compatible version).

    Note: this window contains negative values.
    """

    assert isinstance(n, int) and n >= 1

    coefficients = [1.0 / 4.6402, -1.93 / 4.6402, 1.29 / 4.6402, -0.388 / 4.6402, 0.0322 / 4.6402]

    return cosine_window(n, coefficients, symmetry_flag)


##############################
#                            #
#     TRIANGULAR WINDOWS     #
#                            #
##############################


def triang(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Triangular window.

    triang(1) == [              1.0              ]
    triang(2) == [            0.5 0.5            ]
    triang(3) == [          0.5 1.0 0.5          ]
    triang(4) == [      0.25 0.75 0.75 0.25      ]
    triang(5) == [    0.33 0.66 1.0 0.66 0.33    ]
    triang(6) == [ 0.16 0.50 0.83 0.83 0.50 0.16 ]
    """

    assert isinstance(n, int) and n >= 1

    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return triang(n + 1, True)[:-1].copy()

    if n % 2 == 0:
        # Even length:
        # Center values are (1 - 1 / n); extrema are (1 / n).
        return 1 - np.abs(2 * np.arange(n) - (n - 1)) / (n    )
    else:
        # Odd length:
        # Center value is 1; extrema are 2 / (n + 1).
        return 1 - np.abs(2 * np.arange(n) - (n - 1)) / (n + 1)


def bartlett(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Bartlett window.

    Note: The first and last values of the Bartlett window are zero, for n >= 2.

    bartlett(1) == [           1.0           ]
    bartlett(2) == [         0.0 0.0         ]
    bartlett(3) == [       0.0 1.0 0.0       ]
    bartlett(4) == [    0.0 0.66 0.66 0.0    ]
    bartlett(5) == [   0.0 0.5 1.0 0.5 0.0   ]
    bartlett(6) == [ 0.0 0.4 0.8 0.8 0.4 0.0 ]
    """

    assert isinstance(n, int) and n >= 1

    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return bartlett(n + 1, True)[:-1].copy()

    # Center value is 1 for odd length, 1 - 1 / (n - 1) for even length.
    # Extrema are 0.

    return 1 - np.abs(2 * np.arange(n) - (n - 1)) / (n - 1)

#########################
#                       #
#     OTHER WINDOWS     #
#                       #
#########################

def cosine_scipy(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Cosine window, as defined in scipy."""

    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return cosine_scipy(n + 1, True)[:-1].copy()

    return np.cos((np.arange(n) - (n - 1) / 2) / n * np.pi)


def expwin(n: int, a: float, canonical_flag: bool=False) -> np.ndarray:
    """The exponential window, as implemented in Octave."""

    assert isinstance(n, int) and n >= 1

    if n < 3:
        return np.ones(n)

    # If given SLL then estimate alpha
    if a <= -50.0:
        a = -1.085e-4*a*a - 0.1506*a - 0.304
    elif a <= -13.26:
        a = -1.513e-3*a*a - 0.2809*a - 3.398
    elif a < 0.0:
        a = 0.0

    if canonical_flag:
        T = np.arange(1-n, 0, 2) / (n - 1)
    else:
        T = np.arange(1-n, 0, 2) / (n - 0)

    T = np.sqrt(1.0 - T**2) - 1.0

    W = np.exp(a*T)

    if n % 2 == 0:
        W = np.concatenate((W, W[::-1]))
    else:
        W = np.concatenate((W, [1.0], W[::-1]))

    return W


def poisswin(n: int, a: float) -> np.ndarray:
    """The Poisson window, as implemented in Octave.

    Similar to the expwin window.
    """

    assert isinstance(n, int) and n >= 1

    if n < 3:
        return np.ones(n)

    T = (np.arange((n + 1) // 2) * 2 - (n - 1) // 2 * 2) / (n - 1)
    W = np.exp(a*T)

    if n % 2 == 0:
        W = np.concat((W, W[::-1]))
    else:
        W = np.concat((W, W[-2::-1]))

    return W

def barthannwin(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Modified Bartlett-Hann window."""

    assert isinstance(n, int) and n >= 1

    if n == 1:
        return np.ones(n)

    if not symmetry_flag:
        return barthannwin(n + 1, True)[:-1].copy()

    x = np.abs(np.arange(n) / (n - 1) - 0.5)

    return 0.62 - 0.48 * x + 0.38 * np.cos(2 * np.pi * x)


def bohmanwin(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Bohmann window."""

    assert isinstance(n, int) and n >= 1

    if n == 1:
        return np.ones(n)

    if not symmetry_flag:
        return bohmanwin(n + 1, True)[:-1].copy()

    x = np.abs(2 * np.arange(n) - (n - 1)) / (n - 1)

    return (1 - x) * np.cos(np.pi * x) + np.sin(np.pi * x) / np.pi


def gaussian(n: int, alpha: float = 1.0):
    """Gaussian window, for reproducing the gauss() function in Octave.

    Similar to the gausswin()
    """

    assert isinstance(n, int) and n >= 1

    # x varies from -(n-1)/2 .. +(n-1)/2
    x = np.arange(n)-(n-1)/2

    return np.exp(-0.5*(alpha * x)**2)


def gaussian_scipy(n: int, std: float, symmetry_flag: bool=True):
    """Gaussian window, for reproducing the scipy.signal.windows.gaussian() function from SciPy."""

    assert isinstance(n, int) and n >= 1

    if n == 1:
        return np.ones(n)

    if not symmetry_flag:
        return gaussian_scipy(n + 1, std, True)[:-1].copy()

    # x varies from -(n-1)/2 .. +(n-1)/2
    x = np.arange(n)-(n-1)/2

    return np.exp(-x**2 / (2.0 * std * std))


def gausswin(n: int, alpha: float = 2.5) -> np.ndarray:
    """Gaussian window.

    The parameter for the gausswin() function is different for the Matlab/Octave,
    and SciPy versions of this function:

    - Matlab and Octave use "alpha"
    - Scipy uses "std".

    Matlab vs SciPy: alpha * std == (n - 1) / 2

    In this implementation, we follow the Matlab/Octave convention.
    """

    assert isinstance(n, int) and n >= 1

    # Special case for n == 1, otherwise we'd divide by zero.
    if n == 1:
        return np.ones(1)

    # x varies from -1 .. +1.
    x = (2 * np.arange(n) - (n - 1)) / (n - 1)

    return np.exp( -0.5 * (alpha * x) ** 2)


def parzenwin(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Parzen window.

    This is an approximation of the Gaussian window.

    The Gaussian shape is approximated by two different polynomials,
    one for x < 0.5 and one for x > 0.5. At x == 0.5, the polynomials meet.

    The minimum value of the two polynomials is taken.
    """

    assert isinstance(n, int) and n >= 1

    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return parzenwin(n + 1, True)[:-1].copy()

    x = np.abs(2 * np.arange(n) - (n - 1)) / n

    return np.minimum(1 - 6 * x * x + 6 * x * x * x, 2 * (1 - x) ** 3)


def lanczos(n: int, symmetry_flag: bool=True) -> np.ndarray:
    """Lanczos or sinc window."""

    # Special case for n == 1, otherwise we'd divide by zero.
    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return lanczos(n + 1, True)[:-1].copy()

    x = (np.arange(n) * 2 - (n - 1)) / (n - 1)
    return np.sinc(x)

def tukeywin(n: int, r: float = 0.5, symmetry_flag: bool=True) -> np.ndarray:
    """Tukey window.

    This window uses a cosine-shaped ramp-up and ramp-down, with an all-one part in the middle.
    The parameter 'r' defines the fraction of the window covered by the ramp-up and ramp-down.

    r <= 0 is identical to a rectangular window.
    r >= 1 is identical to a Hann window.
    """

    assert isinstance(n, int) and n >= 1

    # Special case for n == 1, otherwise we'd divide by zero.
    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return tukeywin(n + 1, r, True)[:-1].copy()

    r = np.clip(r, 0, 1)

    # Prevent division by zero for r == 0.
    if r == 0:
        return np.ones(n)

    return (np.cos(np.maximum(np.abs(np.arange(n) - (n - 1) / 2) * (2 / (n - 1) / r)  - (1 / r - 1), 0) * np.pi) + 1) / 2


def taylorwin(n: int, nbar: int = 4, sll: float = -30.0) -> np.ndarray:
    """Taylor window."""

    assert isinstance(n, int) and n >= 1

    # sll is in dB(power).
    # Calculate the amplification factor, e.g. sll = -60 --> amplification = 1000.0

    # Amplification will be >= 1.0.
    amplification = 10 ** (abs(sll) / 20.0)

    # Which is convenient, since arccosh(x) is only real-valued for (x >= 1).
    a = np.arccosh(amplification) / np.pi

    a2 = a * a

    # Taylor pulse widening (dilation) factor.
    sp2 = (nbar ** 2) / (a2 + (nbar - 0.5) ** 2)

    w = np.ones(n)

    x = (np.arange(n) + 0.5) / n

    for m in range(1, nbar):

        # Calculate Fm as a function of: m, sp2, a

        numerator = 1.0
        denominator = 1.0
        for i in range(1, nbar):
            numerator       *= (1.0 - m ** 2 / (sp2 * (a2 + (i - 0.5) ** 2)))
            if i != m:
                denominator *= (1.0 - m ** 2 / i ** 2)

        Fm = -(numerator / denominator)

        w += Fm * np.cos(2 * np.pi * m * x)

    return w


def kaiser(n: int, beta: float = 0.5, symmetry_flag: bool=True) -> np.ndarray:
    """Kaiser window."""

    assert isinstance(n, int) and n >= 1

    # Special case for n == 1, otherwise we'd divide by zero.
    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return kaiser(n + 1, beta, True)[:-1].copy()

    x = (2 * np.arange(n) - (n - 1)) / (n - 1)

    w = (bessel_i0(beta * np.sqrt(1 - x * x)) / bessel_i0(beta))
    return w


def dft_direct(x: np.ndarray) -> np.ndarray:
    """Direct (slow!) implementation of the discrete Fourier transform."""
    n = len(x)
    fx = np.zeros(n, np.complex128)
    for i in range(n):
        fxi = 0
        for j in range(n):
            fxi += np.exp(-float(i * j) / n * 1.j * 2 * np.pi) * x[j]
        fx[i] = fxi
    return fx


def dft(x: np.ndarray) -> np.ndarray:
    """Discrete Fourier transform (defer to numpy)."""
    return np.fft.fft(x)


def chebwin(n: int, r: float = 100.0, symmetry_flag: bool=True) -> np.ndarray:
    """Chebyshev window."""

    assert isinstance(n, int) and n >= 1

    # Special case for n == 1, otherwise we'd divide by zero.
    if n == 1:
        return np.ones(1)

    if not symmetry_flag:
        return chebwin(n + 1, r, True)[:-1].copy()

    order = n - 1

    # r is in dB(power).
    # Calculate the amplification factor, e.g. r = -60 --> amplification = 1000.0

    amplification = 10.0 ** (abs(r) / 20.0)

    beta = np.cosh(np.arccosh(amplification) / order)

    # Find the window's DFT coefficients
    p = np.zeros(n, dtype = np.complex128)

    # Appropriate IDFT and filling up, depending on even/odd length.

    if n % 2 != 0:

        # Odd length window

        for i in range(n):

            x = beta * np.cos(np.pi * i / n)

            if x > 1:
                p[i] = np.cosh(order * np.arccosh( x))
            elif x < -1:
                p[i] = np.cosh(order * np.arccosh(-x))
            else:
                p[i] = np.cos (order * np.arccos ( x))

        p = np.real(dft(p))

        # Example: n = 11
        #
        # w[0] w[1] w[2] w[3] w[4] w[5] w[6] w[7] w[8] w[9] w[10]
        #
        #                            =
        #
        # p[5] p[4] p[3] p[2] p[1] p[0] p[1] p[2] p[3] p[4] p[5]

        nn = (n + 1) // 2

        w = np.concatenate((p[nn - 1:0:-1], p[0:nn]))

    else:

        # Even length window

        for i in range(n):

            x = beta * np.cos(np.pi * i / n)

            z = np.exp(1.j * np.pi * i / n)

            if x > 1:
                p[i] =    z * np.cosh(order * np.arccosh( x))
            elif x < -1:
                p[i] =   -z * np.cosh(order * np.arccosh(-x))
            else:
                p[i] =    z * np.cos (order * np.arccos ( x))

        p = np.real(dft(p))

        # Example: n = 10
        #
        # w[0] w[1] w[2] w[3] w[4] w[5] w[6] w[7] w[8] w[9]
        #
        #                         =
        #
        # p[5] p[4] p[3] p[2] p[1] p[1] p[2] p[3] p[4] p[5]

        nn = n // 2 + 1
        w = np.concatenate((p[nn - 1:0:-1], p[1:nn]))

    # Normalize the window so the maximum value is 1.
    w /= np.amax(w)

    return w
