
import numpy as np
from enum import Enum

# These functions aim to reproduce the behavior of the Matlab functions with the same name.

# The implementations below use only basic NumPy functions; the aim is to have reference implementations
# from which we can easily derive a C implementation.

class SFlag:
    SYMMETRIC = 1
    PERIODIC  = 2

###################################
#                                 #
#     START OF COSINE WINDOWS     #
#                                 #
###################################

def cosine_window(N, coeff, sflag):
    """Generalized cosine window.

        Many window functions described in signal processing literature
        can be written as linear combinations of cosines over the window length.

        Let 'x' be values going from 0 for the first element, and 2*pi for the last element.

        The window can then be written as:

        w = c0 * cos(0 * x) + c1 * cos(1 * x) + c2 * cos(2 * x) + c3 * cos(3 * x) + ...

        Examples of cosine windows implemented below:

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
        The signal processing literature in fact describes many different 'flattop' windows.

        Note 1 : Octave defines the flattopwin coefficients differently, see implementation below.

                 The coefficient values used correspond to:

                 [0.21550795224343777, -0.4159303478298349, 0.2780052583940347, -0.08361708547045386, 0.006939356062238697]

        Note 2 : Octave defines the nuttallwin coefficients differently, see implementation below:

                 The coefficient values used are:

                 [0.355768, -0.487396, 0.144232, -0.012604]
    """

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    # Special case for N == 1, otherwise we'd divide by zero.
    if N <= 1:
        return np.ones(N)

    if sflag == SFlag.SYMMETRIC:
        # The normal, symmetric form of the window.
        x = np.arange(N) * (2 * np.pi / (N - 1))
    else:
        # The periodic for. Calculate the window as if it's one element longer, and discard the last element.
        x = np.arange(N) * (2 * np.pi / (N    ))

    w = np.zeros(N)
    for i in range(len(coeff)):
        w += coeff[i] * np.cos(i * x)

    return w

def rectwin(L):
    """Rectangular window"""

    # Note: this window can also be seen as the simplest of the cosine windows:
    #
    #    return cosine_window(L, [1.0], SFlag.SYMMETRIC)
    #
    # Of course, the implementation given below is faster.

    return np.ones(L) # All values are 1.

def hann(L, sflag = SFlag.SYMMETRIC):
    """Hann (Hanning) window"""

    # Extrema are 0.
    # Center value is 1 for odd length,
    #     0.5 - 0.5 * cos(pi * L / (L - 1)) for even length.

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [0.5, -0.5]

    return cosine_window(L, coeff, sflag)

def hamming(L, sflag = SFlag.SYMMETRIC):
    """Hamming window"""

    # Extrema are 0.08.
    # Center value is 1 for odd length;
    #   0.54 - 0.46 * cos(pi * L / (L - 1)) for even length.

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [0.54, -0.46]

    return cosine_window(L, coeff, sflag)

def blackman(N, sflag = SFlag.SYMMETRIC):
    """Blackman window"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [0.42, -0.5, 0.08]

    return cosine_window(N, [0.42, -0.5, 0.08], sflag)

def blackmanharris(N, sflag = SFlag.SYMMETRIC):
    """Blackman-Harris window"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [0.35875, -0.48829, 0.14128, -0.01168]

    return cosine_window(N, coeff, sflag)

def nuttallwin(N, sflag = SFlag.SYMMETRIC):
    """Nuttall window"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [0.3635819, -0.4891775, 0.1365995, -0.0106411]

    return cosine_window(N, coeff, sflag)

def nuttallwin_octave(N, sflag = SFlag.SYMMETRIC):
    """Nuttall window (Octave version)"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [0.355768, -0.487396, 0.144232, -0.012604]

    return cosine_window(N, coeff, sflag)

def flattopwin(L, sflag = SFlag.SYMMETRIC):
    """Flattop window"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [0.21557895, -0.41663158, 0.277263158, -0.083578947, 0.006947368]

    return cosine_window(L, coeff, sflag)

def flattopwin_octave(L, sflag = SFlag.SYMMETRIC):
    """Flattop window (Octave version)"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    coeff = [1 / 4.6402, -1.93 / 4.6402, 1.29 / 4.6402, -0.388 / 4.6402, 0.0322 / 4.6402]

    return cosine_window(L, coeff, sflag)

##############################
#                            #
#     TRIANGULAR WINDOWS     #
#                            #
##############################

def triang(L):
    """ Triangular window"""
    if L % 2 == 0:
        # Even length
        # Center values are (1 - 1 / L); extrema are (1 / L).
        return 1 - np.abs(2 * np.arange(L) - (L - 1)) / (L    )
    else:
        # Odd length
        # Center value is 1; extrema are 2 / (L + 1)
        return 1 - np.abs(2 * np.arange(L) - (L - 1)) / (L + 1)

def bartlett(L):
    """Bartlett window"""
    # Center value is 1 for odd length, 1 - 1 / (L - 1) for even length.
    # Extrema are 0.

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    return 1 - np.abs(2 * np.arange(L) - (L - 1)) / (L - 1)

#########################
#                       #
#     OTHER WINDOWS     #
#                       #
#########################

def barthannwin(L):
    """Modified Bartlett-Hann window."""

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    x = np.abs(np.arange(L) / (L - 1) - 0.5)

    return 0.62 - 0.48 * x + 0.38 * np.cos(2 * np.pi * x)

def bohmanwin(L):
    """Bohmann window."""

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    x = np.abs(2 * np.arange(L) - (L - 1)) / (L - 1)

    return (1 - x) * np.cos(np.pi * x) + np.sin(np.pi * x) / np.pi

def gausswin(N, Alpha = 2.5):
    """Gaussian window."""
    # Special case for N == 1, otherwise we'd divide by zero.
    if N <= 1:
        return np.ones(N)

    x = np.abs(2 * np.arange(N) - (N - 1)) / (N - 1)

    return np.exp( -0.5 * (Alpha * x) ** 2)

def parzenwin(L):
    """The Parzen window.

       This is an approximation of the Gaussian window.
       The Gaussian shape is approximated by two different polynomials, one for x < 0.5 and one for x > 0.5.
       At x == 0.5, the polynomials meet. The minimum value of the two polynomials is taken.
    """

    x = np.abs(2 * np.arange(L) - (L - 1)) / L

    return np.minimum(1 - 6 * x * x + 6 * x * x * x, 2 * (1 - x) ** 3)

def tukeywin(L, r = 0.5):
    """Tukey window.

       This window uses a cosine-shaped ramp-up and ramp-down, with an all-one part in the middle.
       The parameter 'r' defines the fraction of the window covered by the ramp-up and ramp-down.

       r <= 0 is identical to a rectangular window.
       r >= 1 is identical to a Hann window.
    """

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    r = np.clip(r, 0, 1)

    # Prevent division by zero for r == 0.
    if r == 0:
        return np.ones(L)

    return (np.cos(np.maximum(np.abs(np.arange(L) - (L - 1) / 2) * (2 / (L - 1) / r)  - (1 / r - 1), 0) * np.pi) + 1) / 2

def taylorwin(n, nbar = 4, sll = -30.0):
    """Taylor window."""

    a = np.arccosh(10 ** (-sll / 20.0)) / np.pi

    # Taylor pulse widening (dilation) factor.
    sp2 = (nbar ** 2) / (a ** 2 + (nbar - 0.5) ** 2)

    summation = 0

    xi = (np.arange(n) - ((n - 1) / 2))  / n

    for m in range(1, nbar):

        Num = 1.0
        for i in range(1, nbar):
            Num *= (1 - m ** 2 / sp2 / (a ** 2 + (i - 0.5) ** 2))

        Den = 1.0
        for i in range(1, nbar):
            if i != m:
                Den *= (1 - m**2 / i**2)

        Fm = ((-1) ** (m + 1) * Num) / (2 * Den)

        summation = Fm * np.cos(2 * np.pi * m * xi) + summation

    w = np.ones(n) + 2 * summation

    return w

def bessel_i0(x):
    import scipy.special
    return scipy.special.i0(x)

def kaiser(L, beta = 0.5):
    """Kaiser window."""

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    n = np.arange(0, L)
    alpha = (L - 1) / 2.0
    w = (bessel_i0(beta * np.sqrt(1 - ((n - alpha) / alpha) ** 2.0)) / bessel_i0(beta))
    return w

def chebwin(L, r = 100.0):
    """Chebyshev window."""

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    order = L - 1

    beta = np.cosh(1.0 / order * np.arccosh(10 ** (np.abs(r) / 20.)))

    x = beta * np.cos(np.pi * np.arange(L) / L)

    # Find the window's DFT coefficients
    p = np.zeros(L)
    for i in range(L):
        if x[i] > 1:
            p[i] = np.cosh(order * np.arccosh(x[i]))
        elif x[i] < -1:
            p[i] = (1 - 2 * (order % 2)) * np.cosh(order * np.arccosh(-x[i]))
        else:
            p[i] = np.cos(order * np.arccos(x[i]))

    # Appropriate IDFT and filling up, depending on even/odd length.

    if L % 2 != 0:
        # Odd length
        w = np.real(np.fft.fft(p))
        n = (L + 1) // 2
        w = w[:n]
        w = np.concatenate((w[n - 1:0:-1], w))
    else:
        # Even length
        p = p * np.exp(1.j * np.pi / L * np.arange(L))
        w = np.real(np.fft.fft(p))
        n = L // 2 + 1
        w = np.concatenate((w[n - 1:0:-1], w[1:n]))

    # Normalize window so the maximum value is 1.
    w /= np.amax(w)

    return w
