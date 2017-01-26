
import numpy as np
from enum import Enum

# These functions aim to reproduce the behavior of the Matlab functions with the same name.

# The implementations below use only basic NumPy functions; the aim is to have reference implementations
# from which we can easily derive a C implementation.

class SFlag:
    SYMMETRIC = 1
    PERIODIC  = 2

def rectwin(L):
    """Rectangular window"""

    # All values are 1.

    return np.ones(L)

def hann(L, sflag = SFlag.SYMMETRIC):
    """Hann (Hanning) window"""

    # Extrema are 0.
    # Center value is 1 for odd length,
    #     0.5 - 0.5 * cos(pi * L / (L - 1)) for even length.

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    # Special case for L == 1, otherwise we'd divide by zero.
    # Note that this case also changes the 'periodic' behavior:
    # we return a single 1 value instead of a single value 0.00.

    if L == 1:
        return np.ones(1)

    if sflag == SFlag.PERIODIC:
        return hann(L + 1)[:-1]

    return 0.5 * (1 - np.cos(2 * np.pi * np.arange(L) / (L - 1)))

def hamming(L, sflag = SFlag.SYMMETRIC):
    """Hamming window"""

    # Extrema are 0.08.
    # Center value is 1 for odd length;
    #   0.54 - 0.46 * cos(pi * L / (L - 1)) for even length.

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    # Special case for L == 1, otherwise we'd divide by zero.
    # Note that this case also changes the 'periodic' behavior:
    # we return a single 1 value instead of a single value 0.08.

    if L == 1:
        return np.ones(1)

    if sflag == SFlag.PERIODIC:
        return hamming(L + 1)[:-1]

    return 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(L) / (L - 1))

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
    if L == 1:
        return np.ones(1)

    return 1 - np.abs(2 * np.arange(L) - (L - 1)) / (L - 1)

def barthannwin(L):
    """Modified Bartlett-Hann window"""

    # Special case for L == 1, otherwise we'd divide by zero.
    if L == 1:
        return np.ones(1)

    x = np.abs(np.arange(L) / (L - 1) - 0.5)

    return 0.62 - 0.48 * x + 0.38 * np.cos(2 * np.pi * x)

def nuttallwin(N, sflag = SFlag.SYMMETRIC):
    """Nuttall window"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    # Special case for N == 1, otherwise we'd divide by zero.
    if N == 1:
        return np.ones(1)

    if sflag == SFlag.PERIODIC:
        return nuttallwin(N + 1)[:-1]

    a = [0.3635819, 0.4891775, 0.1365995, 0.0106411]
    n = np.arange(0, N)
    fac = n * 2 * np.pi / (N - 1.0)
    w = (a[0] - a[1] * np.cos(fac) + a[2] * np.cos(2 * fac) - a[3] * np.cos(3 * fac))
    return w

def blackman(N, sflag = SFlag.SYMMETRIC):
    """Blackman window"""
    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    # Special case for N == 1, otherwise we'd divide by zero.
    if N == 1:
        return np.ones(1)

    if sflag == SFlag.PERIODIC:
        return blackman(N + 1)[:-1]

    x = np.arange(N) * (2 * np.pi / (N - 1))

    return 0.42 - 0.5 * np.cos(x) + 0.08 * np.cos(2.0 * x)

def blackmanharris(N, sflag = SFlag.SYMMETRIC):
    """Blackman-Harris window"""

    assert sflag in [SFlag.SYMMETRIC, SFlag.PERIODIC]

    # Special case for N == 1, otherwise we'd divide by zero.
    if N == 1:
        return np.ones(1)

    if sflag == SFlag.PERIODIC:
        return blackmanharris(N + 1)[:-1]

    x = np.arange(N) * (2 * np.pi / (N - 1))

    return 0.35875 - 0.48829 * np.cos(x) + 0.14128 * np.cos(2 * x) - 0.01168 * np.cos(3 * x)

def bohmanwin(L):

    # Special case for L == 1, otherwise we'd divide by zero.
    if L == 1:
        return np.ones(1)

    x = np.abs(2 * np.arange(L) - (L - 1)) / (L - 1)

    return (1 - x) * np.cos(np.pi * x) + np.sin(np.pi * x) / np.pi

def gausswin(N, Alpha = 2.5):

    # Special case for N == 1, otherwise we'd divide by zero.
    if N == 1:
        return np.ones(1)

    x = np.abs(2 * np.arange(N) - (N - 1)) / (N - 1)
    return np.exp( -0.5 * (Alpha * x) ** 2)

def tukeywin(L, r = 0.5):

    # Special case for L == 1, otherwise we'd divide by zero.
    if L == 1:
        return np.ones(1)

    r = np.clip(r, 0, 1)

    # Prevent division by zero for r == 0.
    if r == 0:
        return np.ones(L)

    return (np.cos(np.maximum(np.abs(np.arange(L) - (L - 1) / 2) * (2 / (L - 1) / r)  - (1 / r - 1), 0) * np.pi) + 1) / 2

def parzenwin(L):

    x = np.abs(2 * np.arange(L) - (L - 1)) / L

    return np.minimum(1 - 6 * x * x + 6 * x * x * x, 2 * (1 - x) ** 3)

def flattopwin(L, sflag = True):

    # Special case for L == 1, otherwise we'd divide by zero.
    if L == 1:
        return np.ones(1)

    if sflag == SFlag.PERIODIC:
        return flattopwin(L + 1)[:-1]

    a0 = 0.21557895
    a1 = 0.41663158
    a2 = 0.277263158
    a3 = 0.083578947
    a4 = 0.006947368

    a = [a0, a1, a2, a3, a4]

    n = np.arange(0, L)
    fac = n * 2 * np.pi / (L - 1.0)
    w = (a[0] - a[1] * np.cos(fac) + a[2] * np.cos(2 * fac) - a[3] * np.cos(3 * fac) + a[4] * np.cos(4 * fac))

    return w


def taylorwin(n, nbar = 4, sll = -30.0):

    a = np.arccosh(10 ** (-sll / 20.0)) / np.pi

    # Taylor pulse widening (dilation) factor.
    sp2 = (nbar ** 2) / (a ** 2 + (nbar - 0.5) ** 2)

    summation = 0

    xi = (np.arange(n) - ((n - 1) / 2))  / n

    for m in range(1, nbar):

        Num = 1.0
        for i in range(1, nbar):
            Num *= 1 - m ** 2 / sp2 / (a ** 2 + (i - 0.5) ** 2)

        Den = 1.0
        for i in range(1, nbar):
            if i != m:
                Den *= 1 - m**2 / i**2

        Fm = ((-1) ** (m+1) *Num) / (2 * Den)

        summation = Fm * np.cos(2 * np.pi * m * xi) + summation

    w = np.ones(n) + 2 * summation

    return w
