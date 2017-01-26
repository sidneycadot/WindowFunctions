#! /usr/bin/env python3

import scipy.signal
import numpy as np
from collections import OrderedDict
import window_functions

class ReferenceWindowFunction:

    def __init__(self):
        self.windows = OrderedDict()

    def set_value(self, M, i, value):

        if M in self.windows:
            window = self.windows[M]
        else:
            window = np.zeros(M) + np.nan # initialize window with NaN values.
            self.windows[M] = window

        i -= 1 # Matlab/octave use 1-based indexing, but we use 0-based indexing.

        assert (0 <= i < M)

        # Check that this window is new.
        assert np.isnan(window[i])
        window[i] = value

def read_reference_file(filename):

    ref_winfuncs = OrderedDict()

    for line in open(filename):

        (name, M, i, value) = line.split()
        M = int(M)
        i = int(i)
        value = float(value)

        if name in ref_winfuncs:
            ref_winfunc = ref_winfuncs[name]
        else:
            ref_winfunc = ReferenceWindowFunction()
            ref_winfuncs[name] = ref_winfunc

        ref_winfunc.set_value(M, i, value)

    # Verify that all reference winfunc values are finite (not NaN).

    for ref_winfunc in ref_winfuncs.values():
        for w in ref_winfunc.windows.values():
            assert np.isfinite(w).all()

    # Return all reference window functions
    return ref_winfuncs

def make_winfunc_map():

    # Python equivalents for the window types we find in the reference files ...

    winfunc_map = OrderedDict()

    # DONE:
    #
    #   barthannwin
    #   bartlett
    #   blackman
    #   blackmanharris
    #   bohmanwin
    #   gausswin
    #   hamming
    #   hann
    #   nuttallwin
    #   parzenwin
    #   rectwin
    #   triang
    #   tukeywin
    #
    # TODO:
    #
    #    chebwin    (SciPy version works.)
    #    kaiser     (SciPy version works.)
    #    taylor

    winfunc_map[ "barthannwin"              ] = lambda M : window_functions.barthannwin(M)
    winfunc_map[ "bartlett"                 ] = lambda M : window_functions.bartlett(M)
    winfunc_map[ "blackman"                 ] = lambda M : window_functions.blackman(M)
    winfunc_map[ "blackman_periodic"        ] = lambda M : window_functions.blackman(M, window_functions.SFlag.PERIODIC)
    winfunc_map[ "blackman_symmetric"       ] = lambda M : window_functions.blackman(M, window_functions.SFlag.SYMMETRIC)
    winfunc_map[ "blackmanharris"           ] = lambda M : window_functions.blackmanharris(M)
    winfunc_map[ "blackmanharris_periodic"  ] = lambda M : window_functions.blackmanharris(M, window_functions.SFlag.PERIODIC)
    winfunc_map[ "blackmanharris_symmetric" ] = lambda M : window_functions.blackmanharris(M, window_functions.SFlag.SYMMETRIC)
    winfunc_map[ "bohmanwin"                ] = lambda M : window_functions.bohmanwin(M)
    winfunc_map[ "chebwin"                  ] = lambda M : scipy.signal.chebwin(M, 100.0)
    winfunc_map[ "chebwin_100p0"            ] = lambda M : scipy.signal.chebwin(M, 100.0)
    winfunc_map[ "chebwin_120p0"            ] = lambda M : scipy.signal.chebwin(M, 120.0)
    winfunc_map[ "flattopwin"               ] = lambda M : window_functions.flattopwin(M)
    winfunc_map[ "flattopwin_periodic"      ] = lambda M : window_functions.flattopwin(M, window_functions.SFlag.PERIODIC)
    winfunc_map[ "flattopwin_symmetric"     ] = lambda M : window_functions.flattopwin(M, window_functions.SFlag.SYMMETRIC)
    winfunc_map[ "gausswin"                 ] = lambda M : window_functions.gausswin(M)
    winfunc_map[ "gausswin_2p5"             ] = lambda M : window_functions.gausswin(M, 2.5)
    winfunc_map[ "gausswin_3p2"             ] = lambda M : window_functions.gausswin(M, 3.2)
    winfunc_map[ "hamming"                  ] = lambda M : window_functions.hamming(M)
    winfunc_map[ "hamming_periodic"         ] = lambda M : window_functions.hamming(M, window_functions.SFlag.PERIODIC)
    winfunc_map[ "hamming_symmetric"        ] = lambda M : window_functions.hamming(M, window_functions.SFlag.SYMMETRIC)
    winfunc_map[ "hann"                     ] = lambda M : window_functions.hann(M)
    winfunc_map[ "hann_periodic"            ] = lambda M : window_functions.hann(M, window_functions.SFlag.PERIODIC)
    winfunc_map[ "hann_symmetric"           ] = lambda M : window_functions.hann(M, window_functions.SFlag.SYMMETRIC)
    winfunc_map[ "kaiser"                   ] = lambda M : scipy.signal.kaiser(M, 0.5)
    winfunc_map[ "kaiser_0p5"               ] = lambda M : scipy.signal.kaiser(M, 0.5)
    winfunc_map[ "kaiser_0p8"               ] = lambda M : scipy.signal.kaiser(M, 0.8)
    winfunc_map[ "nuttallwin"               ] = lambda M : window_functions.nuttallwin(M)
    winfunc_map[ "nuttallwin_periodic"      ] = lambda M : window_functions.nuttallwin(M, window_functions.SFlag.PERIODIC)
    winfunc_map[ "nuttallwin_symmetric"     ] = lambda M : window_functions.nuttallwin(M, window_functions.SFlag.SYMMETRIC)
    winfunc_map[ "parzenwin"                ] = lambda M : window_functions.parzenwin(M)
    winfunc_map[ "rectwin"                  ] = lambda M : window_functions.rectwin(M)
    winfunc_map[ "taylorwin"                ] = lambda M : window_functions.taylorwin(M)
    winfunc_map[ "taylorwin_4"              ] = lambda M : window_functions.taylorwin(M, 4)
    winfunc_map[ "taylorwin_5"              ] = lambda M : window_functions.taylorwin(M, 5)
    winfunc_map[ "taylorwin_6"              ] = lambda M : window_functions.taylorwin(M, 6)
    winfunc_map[ "taylorwin_4_m20"          ] = lambda M : window_functions.taylorwin(M, 4, -20)
    winfunc_map[ "taylorwin_5_m20"          ] = lambda M : window_functions.taylorwin(M, 5, -20)
    winfunc_map[ "taylorwin_6_m20"          ] = lambda M : window_functions.taylorwin(M, 6, -20)
    winfunc_map[ "taylorwin_4_m30"          ] = lambda M : window_functions.taylorwin(M, 4, -30)
    winfunc_map[ "taylorwin_5_m30"          ] = lambda M : window_functions.taylorwin(M, 5, -30)
    winfunc_map[ "taylorwin_6_m30"          ] = lambda M : window_functions.taylorwin(M, 6, -30)
    winfunc_map[ "taylorwin_4_m40"          ] = lambda M : window_functions.taylorwin(M, 4, -40)
    winfunc_map[ "taylorwin_5_m40"          ] = lambda M : window_functions.taylorwin(M, 5, -40)
    winfunc_map[ "taylorwin_6_m40"          ] = lambda M : window_functions.taylorwin(M, 6, -40)
    winfunc_map[ "triang"                   ] = lambda M : window_functions.triang(M)
    winfunc_map[ "tukeywin"                 ] = lambda M : window_functions.tukeywin(M)
    winfunc_map[ "tukeywin_0p0"             ] = lambda M : window_functions.tukeywin(M, 0.0)
    winfunc_map[ "tukeywin_0p2"             ] = lambda M : window_functions.tukeywin(M, 0.2)
    winfunc_map[ "tukeywin_0p5"             ] = lambda M : window_functions.tukeywin(M, 0.5)
    winfunc_map[ "tukeywin_0p8"             ] = lambda M : window_functions.tukeywin(M, 0.8)
    winfunc_map[ "tukeywin_1p0"             ] = lambda M : window_functions.tukeywin(M, 1.0)

    return winfunc_map

def check_reference_window_functions(ref_winfuncs, winfunc_map):

    for (ref_winfunc_name, ref_winfunc) in ref_winfuncs.items():
        if ref_winfunc_name in winfunc_map:
            python_func = winfunc_map[ref_winfunc_name]
            ok = True
            worsterr = 0.0
            for (M, ref_values) in ref_winfunc.windows.items():
                python_values = python_func(M)
                err = np.abs(ref_values - python_values)
                maxerr = max(err) # Maximum error for this window size
                if maxerr > 1e-12:
                    print("Bad window value: {} M = {} maxerr = {}".format(ref_winfunc_name, M, maxerr))
                    print("ref_values:", ref_values)
                    print("python_values:", python_values)
                    ok = False
                worsterr = max(maxerr, worsterr) # Maximum error across all window sizes.
            if ok:
                print("Perfect correspondence ........ : {:30} (worsterr = {:10.10g})".format(ref_winfunc_name, worsterr))
        else:
            print("Not found in winfunc_map ...... : {:30}".format(ref_winfunc_name))
    print()

if True:
    print()
    print("*** CHECK MATLAB REFERENCE VALUES ***")
    print()
    ref_winfuncs = read_reference_file("reference/matlab_windows.txt")

    winfunc_map = make_winfunc_map()
    check_reference_window_functions(ref_winfuncs, winfunc_map)

if False:
    print()
    print("*** CHECK OCTAVE REFERENCE VALUES ***")
    print()
    ref_winfuncs = read_reference_file("reference/octave_windows.txt")

    winfunc_map = make_winfunc_map()
    check_reference_window_functions(ref_winfuncs, winfunc_map)
