#! /usr/bin/env python3

import numpy as np
from collections import OrderedDict
import window_functions as wf

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

def make_matlab_winfunc_map():

    # Python equivalents for the window types we find in the Matplab reference files.

    winfunc_map = OrderedDict()

    winfunc_map[ "barthannwin"              ] = lambda M : wf.barthannwin(M)
    winfunc_map[ "bartlett"                 ] = lambda M : wf.bartlett(M)
    winfunc_map[ "blackman"                 ] = lambda M : wf.blackman(M)
    winfunc_map[ "blackman_periodic"        ] = lambda M : wf.blackman(M, False)
    winfunc_map[ "blackman_symmetric"       ] = lambda M : wf.blackman(M, True)
    winfunc_map[ "blackmanharris"           ] = lambda M : wf.blackmanharris(M)
    winfunc_map[ "blackmanharris_periodic"  ] = lambda M : wf.blackmanharris(M, False)
    winfunc_map[ "blackmanharris_symmetric" ] = lambda M : wf.blackmanharris(M, True)
    winfunc_map[ "bohmanwin"                ] = lambda M : wf.bohmanwin(M)
    winfunc_map[ "chebwin"                  ] = lambda M : wf.chebwin(M)
    winfunc_map[ "chebwin_100p0"            ] = lambda M : wf.chebwin(M, 100.0)
    winfunc_map[ "chebwin_120p0"            ] = lambda M : wf.chebwin(M, 120.0)
    winfunc_map[ "flattopwin"               ] = lambda M : wf.flattopwin(M)
    winfunc_map[ "flattopwin_periodic"      ] = lambda M : wf.flattopwin(M, False)
    winfunc_map[ "flattopwin_symmetric"     ] = lambda M : wf.flattopwin(M, True)
    winfunc_map[ "gausswin"                 ] = lambda M : wf.gausswin(M)
    winfunc_map[ "gausswin_2p5"             ] = lambda M : wf.gausswin(M, 2.5)
    winfunc_map[ "gausswin_3p2"             ] = lambda M : wf.gausswin(M, 3.2)
    winfunc_map[ "hamming"                  ] = lambda M : wf.hamming(M)
    winfunc_map[ "hamming_periodic"         ] = lambda M : wf.hamming(M, False)
    winfunc_map[ "hamming_symmetric"        ] = lambda M : wf.hamming(M, True)
    winfunc_map[ "hann"                     ] = lambda M : wf.hann(M)
    winfunc_map[ "hann_periodic"            ] = lambda M : wf.hann(M, False)
    winfunc_map[ "hann_symmetric"           ] = lambda M : wf.hann(M, True)
    winfunc_map[ "kaiser"                   ] = lambda M : wf.kaiser(M)
    winfunc_map[ "kaiser_0p5"               ] = lambda M : wf.kaiser(M, 0.5)
    winfunc_map[ "kaiser_0p8"               ] = lambda M : wf.kaiser(M, 0.8)
    winfunc_map[ "nuttallwin"               ] = lambda M : wf.nuttallwin(M)
    winfunc_map[ "nuttallwin_periodic"      ] = lambda M : wf.nuttallwin(M, False)
    winfunc_map[ "nuttallwin_symmetric"     ] = lambda M : wf.nuttallwin(M, True)
    winfunc_map[ "parzenwin"                ] = lambda M : wf.parzenwin(M)
    winfunc_map[ "rectwin"                  ] = lambda M : wf.rectwin(M)
    winfunc_map[ "taylorwin"                ] = lambda M : wf.taylorwin(M)
    winfunc_map[ "taylorwin_4"              ] = lambda M : wf.taylorwin(M, 4)
    winfunc_map[ "taylorwin_5"              ] = lambda M : wf.taylorwin(M, 5)
    winfunc_map[ "taylorwin_6"              ] = lambda M : wf.taylorwin(M, 6)
    winfunc_map[ "taylorwin_4_m20"          ] = lambda M : wf.taylorwin(M, 4, -20)
    winfunc_map[ "taylorwin_5_m20"          ] = lambda M : wf.taylorwin(M, 5, -20)
    winfunc_map[ "taylorwin_6_m20"          ] = lambda M : wf.taylorwin(M, 6, -20)
    winfunc_map[ "taylorwin_4_m30"          ] = lambda M : wf.taylorwin(M, 4, -30)
    winfunc_map[ "taylorwin_5_m30"          ] = lambda M : wf.taylorwin(M, 5, -30)
    winfunc_map[ "taylorwin_6_m30"          ] = lambda M : wf.taylorwin(M, 6, -30)
    winfunc_map[ "taylorwin_4_m40"          ] = lambda M : wf.taylorwin(M, 4, -40)
    winfunc_map[ "taylorwin_5_m40"          ] = lambda M : wf.taylorwin(M, 5, -40)
    winfunc_map[ "taylorwin_6_m40"          ] = lambda M : wf.taylorwin(M, 6, -40)
    winfunc_map[ "triang"                   ] = lambda M : wf.triang(M)
    winfunc_map[ "tukeywin"                 ] = lambda M : wf.tukeywin(M)
    winfunc_map[ "tukeywin_0p0"             ] = lambda M : wf.tukeywin(M, 0.0)
    winfunc_map[ "tukeywin_0p2"             ] = lambda M : wf.tukeywin(M, 0.2)
    winfunc_map[ "tukeywin_0p5"             ] = lambda M : wf.tukeywin(M, 0.5)
    winfunc_map[ "tukeywin_0p8"             ] = lambda M : wf.tukeywin(M, 0.8)
    winfunc_map[ "tukeywin_1p0"             ] = lambda M : wf.tukeywin(M, 1.0)

    return winfunc_map

def make_octave_winfunc_map():

    # Python equivalents for the window types we find in the Octave reference files.

    # Differences with the Matlab functions:
    #
    #   * The flattop window is defined differently (a cosine window with slightly different doefficients).
    #   * The Nutall window is defined differently (a cosine window with slightly different doefficients).
    #   * The Gauss window is defined differently. The parameter is scaled by a factor M/(M - 1) relative to the equivalant matlab function.
    #   * Octave does not implement the Taylow window.

    winfunc_map = OrderedDict()

    winfunc_map[ "barthannwin"              ] = lambda M : wf.barthannwin(M)
    winfunc_map[ "bartlett"                 ] = lambda M : wf.bartlett(M)
    winfunc_map[ "blackman"                 ] = lambda M : wf.blackman(M)
    winfunc_map[ "blackman_periodic"        ] = lambda M : wf.blackman(M, False)
    winfunc_map[ "blackman_symmetric"       ] = lambda M : wf.blackman(M, True)
    winfunc_map[ "blackmanharris"           ] = lambda M : wf.blackmanharris(M)
    winfunc_map[ "blackmanharris_periodic"  ] = lambda M : wf.blackmanharris(M, False)
    winfunc_map[ "blackmanharris_symmetric" ] = lambda M : wf.blackmanharris(M, True)
    winfunc_map[ "bohmanwin"                ] = lambda M : wf.bohmanwin(M)
    winfunc_map[ "chebwin"                  ] = lambda M : wf.chebwin(M)
    winfunc_map[ "chebwin_100p0"            ] = lambda M : wf.chebwin(M, 100.0)
    winfunc_map[ "chebwin_120p0"            ] = lambda M : wf.chebwin(M, 120.0)
    winfunc_map[ "flattopwin"               ] = lambda M : wf.flattopwin_octave(M)
    winfunc_map[ "flattopwin_periodic"      ] = lambda M : wf.flattopwin_octave(M, False)
    winfunc_map[ "flattopwin_symmetric"     ] = lambda M : wf.flattopwin_octave(M, True)
    winfunc_map[ "gausswin"                 ] = lambda M : wf.gausswin(M, 2.5 * (M - 1) / M)
    winfunc_map[ "gausswin_2p5"             ] = lambda M : wf.gausswin(M, 2.5 * (M - 1) / M)
    winfunc_map[ "gausswin_3p2"             ] = lambda M : wf.gausswin(M, 3.2 * (M - 1) / M)
    winfunc_map[ "hamming"                  ] = lambda M : wf.hamming(M)
    winfunc_map[ "hamming_periodic"         ] = lambda M : wf.hamming(M, False)
    winfunc_map[ "hamming_symmetric"        ] = lambda M : wf.hamming(M, True)
    winfunc_map[ "hann"                     ] = lambda M : wf.hann(M)
    winfunc_map[ "hann_periodic"            ] = lambda M : wf.hann(M, False)
    winfunc_map[ "hann_symmetric"           ] = lambda M : wf.hann(M, True)
    winfunc_map[ "kaiser"                   ] = lambda M : wf.kaiser(M)
    winfunc_map[ "kaiser_0p5"               ] = lambda M : wf.kaiser(M, 0.5)
    winfunc_map[ "kaiser_0p8"               ] = lambda M : wf.kaiser(M, 0.8)
    winfunc_map[ "nuttallwin"               ] = lambda M : wf.nuttallwin_octave(M)
    winfunc_map[ "nuttallwin_periodic"      ] = lambda M : wf.nuttallwin_octave(M, False)
    winfunc_map[ "nuttallwin_symmetric"     ] = lambda M : wf.nuttallwin_octave(M, True)
    winfunc_map[ "parzenwin"                ] = lambda M : wf.parzenwin(M)
    winfunc_map[ "rectwin"                  ] = lambda M : wf.rectwin(M)
    winfunc_map[ "triang"                   ] = lambda M : wf.triang(M)
    winfunc_map[ "tukeywin"                 ] = lambda M : wf.tukeywin(M)
    winfunc_map[ "tukeywin_0p0"             ] = lambda M : wf.tukeywin(M, 0.0)
    winfunc_map[ "tukeywin_0p2"             ] = lambda M : wf.tukeywin(M, 0.2)
    winfunc_map[ "tukeywin_0p5"             ] = lambda M : wf.tukeywin(M, 0.5)
    winfunc_map[ "tukeywin_0p8"             ] = lambda M : wf.tukeywin(M, 0.8)
    winfunc_map[ "tukeywin_1p0"             ] = lambda M : wf.tukeywin(M, 1.0)

    return winfunc_map

def check_reference_wf(ref_winfuncs, winfunc_map):

    for (ref_winfunc_name, ref_winfunc) in ref_winfuncs.items():
        if ref_winfunc_name in winfunc_map:
            python_func = winfunc_map[ref_winfunc_name]
            ok = True
            worsterr = 0.0
            for (M, ref_values) in ref_winfunc.windows.items():
                python_values = python_func(M)
                err = np.abs(ref_values - python_values)
                maxerr = max(err) # Maximum error for this window size
                if maxerr > 1e-13:
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

    winfunc_map = make_matlab_winfunc_map()
    check_reference_wf(ref_winfuncs, winfunc_map)

if True:
    print()
    print("*** CHECK OCTAVE REFERENCE VALUES ***")
    print()
    ref_winfuncs = read_reference_file("reference/octave_windows.txt")

    winfunc_map = make_octave_winfunc_map()
    check_reference_wf(ref_winfuncs, winfunc_map)
