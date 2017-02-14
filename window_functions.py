
import numpy as np

# These functions aim to reproduce the behavior of the Matlab functions with the same name.

# The implementations below use only basic NumPy functions; the aim is to have reference implementations
# from which we can easily derive a C implementation.
#
# There are two window types for which this is somewhat tricky:
#
# - The Kaiser window depends on the ability to calculate the modified Bessel function of the first kind I0(x).
# - The Chebyshev window depends on the ability top perform a Discrete Fourier Transform (DFT).

###################################
#                                 #
#     START OF COSINE WINDOWS     #
#                                 #
###################################

def cosine_window(N, coeff, sflag = True):
    """ Generalized cosine window.

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
        The signal processing literature in fact describes many different 'flattop' windows.

        Note 1 : Octave defines the flattopwin coefficients differently, see implementation below.

                 The coefficient values used correspond to:

                 [0.21550795224343777, -0.4159303478298349, 0.2780052583940347, -0.08361708547045386, 0.006939356062238697]

        Note 2 : Octave defines the nuttallwin coefficients differently, see implementation below:

                 The coefficient values used are:

                 [0.355768, -0.487396, 0.144232, -0.012604]
    """

    # Special case for N == 1, otherwise we'd divide by zero.
    if N <= 1:
        return np.ones(N)

    if sflag:
        # The normal, symmetric form of the window.
        x = np.arange(N) * (2 * np.pi / (N - 1))
    else:
        # The periodic form. Calculate the window as if it's one element longer, and discard the last element.
        x = np.arange(N) * (2 * np.pi / (N    ))

    w = np.zeros(N)
    for i in range(len(coeff)):
        w += coeff[i] * np.cos(i * x)

    return w

def rectwin(L):
    """ Rectangular window
    """

    # Note: this window can also be seen as the simplest of the cosine windows:
    #
    #    return cosine_window(L, [1.0], True)
    #
    # Of course, the implementation given below is faster.

    return np.ones(L) # All values are 1.

def hann(L, sflag = True):
    """ Hann window
    """

    # Extrema are 0.
    # Center value is 1 for odd length,
    #     0.5 - 0.5 * cos(pi * L / (L - 1)) for even length.

    coeff = [0.5, -0.5]

    return cosine_window(L, coeff, sflag)

def hamming(L, sflag = True):
    """ Hamming window
    """

    # Note that the Hamming window is raised; its extreme values are 0.08.
    #
    # The center value is 1 for odd length;
    # The venter values are 0.54 - 0.46 * cos(pi * L / (L - 1)) for even length.

    coeff = [0.54, -0.46]

    return cosine_window(L, coeff, sflag)

def blackman(N, sflag = True):
    """ Blackman window
    """

    coeff = [0.42, -0.5, 0.08]

    return cosine_window(N, coeff, sflag)

def blackmanharris(N, sflag = True):
    """ Blackman-Harris window

        Note: very similar to the Nuttall window.
    """

    coeff = [0.35875, -0.48829, 0.14128, -0.01168]

    return cosine_window(N, coeff, sflag)

def nuttallwin(N, sflag = True):
    """ Nuttall window

        Note: very similar to the Blackman-Harris window.
    """

    coeff = [0.3635819, -0.4891775, 0.1365995, -0.0106411]

    return cosine_window(N, coeff, sflag)


    return cosine_window(N, coeff, sflag)

def nuttallwin_octave(N, sflag = True):
    """ Nuttall window (Octave version)
    """

    coeff = [0.355768, -0.487396, 0.144232, -0.012604]

    return cosine_window(N, coeff, sflag)

def flattopwin(L, sflag = True):
    """ Flattop window

        Note: this window contains negative entries!
    """

    coeff = [0.21557895, -0.41663158, 0.277263158, -0.083578947, 0.006947368]

    return cosine_window(L, coeff, sflag)

def flattopwin_octave(L, sflag = True):
    """ Flattop window (Octave version)

        Note: this window contains negative values.
    """

    coeff = [1.0 / 4.6402, -1.93 / 4.6402, 1.29 / 4.6402, -0.388 / 4.6402, 0.0322 / 4.6402]

    return cosine_window(L, coeff, sflag)

##############################
#                            #
#     TRIANGULAR WINDOWS     #
#                            #
##############################

def triang(L):
    """ Triangular window

        triang(1) == [              1.0              ]
        triang(2) == [            0.5 0.5            ]
        triang(3) == [          0.5 1.0 0.5          ]
        triang(4) == [      0.25 0.75 0.75 0.25      ]
        triang(5) == [    0.33 0.66 1.0 0.66 0.33    ]
        triang(6) == [ 0.16 0.50 0.83 0.83 0.50 0.16 ]
    """

    if L % 2 == 0:
        # Even length:
        # Center values are (1 - 1 / L); extrema are (1 / L).
        return 1 - np.abs(2 * np.arange(L) - (L - 1)) / (L    )
    else:
        # Odd length:
        # Center value is 1; extrema are 2 / (L + 1).
        return 1 - np.abs(2 * np.arange(L) - (L - 1)) / (L + 1)

def bartlett(L):
    """ Bartlett window

        bartlett(1) == [           1.0           ]
        bartlett(2) == [         0.0 0.0         ]
        bartlett(3) == [       0.0 1.0 0.0       ]
        bartlett(4) == [    0.0 0.66 0.66 0.0    ]
        bartlett(5) == [   0.0 0.5 1.0 0.5 0.0   ]
        bartlett(6) == [ 0.0 0.4 0.8 0.8 0.4 0.0 ]
    """

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
    """ Modified Bartlett-Hann window.
    """

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    x = np.abs(np.arange(L) / (L - 1) - 0.5)

    return 0.62 - 0.48 * x + 0.38 * np.cos(2 * np.pi * x)

def bohmanwin(L):
    """Bohmann window.
    """

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    x = np.abs(2 * np.arange(L) - (L - 1)) / (L - 1)

    return (1 - x) * np.cos(np.pi * x) + np.sin(np.pi * x) / np.pi

def gausswin(N, Alpha = 2.5):
    """ Gaussian window.

        The parameter for the gausswin() function is different for the Matlab, Octave,
        and SciPy versions of this function:

        - Matlab uses "Alpha";
        - Octave uses "A";
        - Scipy uses "std".

        Matlab vs SciPy:     Alpha * std == (N - 1) / 2

        Matlab vs Octave:    Alpha * N == A * (N - 1)

        In this implementation, we follow the Matlab convention.
    """

    # Special case for N == 1, otherwise we'd divide by zero.
    if N <= 1:
        return np.ones(N)

    x = np.abs(2 * np.arange(N) - (N - 1)) / (N - 1)

    return np.exp( -0.5 * (Alpha * x) ** 2)

def parzenwin(L):
    """ Parzen window.

        This is an approximation of the Gaussian window.

        The Gaussian shape is approximated by two different polynomials,
        one for x < 0.5 and one for x > 0.5. At x == 0.5, the polynomials meet.

        The minimum value of the two polynomials is taken.
    """

    x = np.abs(2 * np.arange(L) - (L - 1)) / L

    return np.minimum(1 - 6 * x * x + 6 * x * x * x, 2 * (1 - x) ** 3)

def tukeywin(L, r = 0.5):
    """ Tukey window.

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
    """ Taylor window.
    """

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

def kaiser(L, beta = 0.5):
    """ Kaiser window.
    """

    from bessel_i0 import bessel_i0

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    x = (2 * np.arange(L) - (L - 1)) / (L - 1)

    w = (bessel_i0(beta * np.sqrt(1 - x * x)) / bessel_i0(beta))
    return w

def dft_direct(x):
    n = len(x)
    fx = np.zeros(n, np.complex128)
    for i in range(n):
        fxi = 0
        for j in range(n):
            fxi += np.exp(-float(i * j) / n * 1.j * 2 * np.pi) * x[j]
        fx[i] = fxi
    return fx

def dft(x):
    return np.fft.fft(x)

def chebwin(L, r = 100.0):
    """ Chebyshev window.
    """

    # Special case for L == 1, otherwise we'd divide by zero.
    if L <= 1:
        return np.ones(L)

    order = L - 1

    # sll is in dB(power).
    # Calculate the amplification factor, e.g. sll = -60 --> amplification = 1000.0

    amplification = 10.0 ** (abs(r) / 20.0)

    beta = np.cosh(np.arccosh(amplification) / order)

    # Find the window's DFT coefficients
    p = np.zeros(L, dtype = np.complex128)

    # Appropriate IDFT and filling up, depending on even/odd length.

    if L % 2 != 0:

        # Odd length window

        for i in range(L):

            x =  beta * np.cos(np.pi * i / L)

            if x > 1:
                p[i] =    np.cosh(order * np.arccosh( x))
            elif x < -1:
                p[i] =    np.cosh(order * np.arccosh(-x))
            else:
                p[i] =    np.cos (order * np.arccos ( x))

        w = np.real(dft(p))

        # w[0] w[1] w[2] w[3] w[4] w[5] w[6] w[7] w[8] w[9] w[10]
        #
        # w[5] w[4] w[3] w[2] w[1] w[0] w[1] w[2] w[3] w[4] w[5]

        n = (L + 1) // 2

        w = np.concatenate((w[n - 1:0:-1], w[0:n]))

    else:

        # Even length window

        for i in range(L):

            x =  beta * np.cos(np.pi * i / L)

            #z = 1 + 0.j
            z = np.exp(1.j * np.pi * i / L)

            if x > 1:
                p[i] =    z * np.cosh(order * np.arccosh( x))
            elif x < -1:
                p[i] =   -z * np.cosh(order * np.arccosh(-x))
            else:
                p[i] =    z * np.cos (order * np.arccos ( x))

        w = np.real(dft(p))

        # w[0] w[1] w[2] w[3] w[4] w[5] w[6] w[7] w[8] w[9]
        #
        # w[5] w[4] w[3] w[2] w[1] w[1] w[2] w[3] w[4] w[5]

        n = L // 2 + 1
        w = np.concatenate((w[n - 1:0:-1], w[1:n]))

    # Normalize window so the maximum value is 1.
    w /= np.amax(w)

    return w
