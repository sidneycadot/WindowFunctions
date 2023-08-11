#! /usr/bin/env python3

import sys
import contextlib
from typing import Callable

import numpy
import scipy

def dump_window_function(source: str, name: str, f: Callable[[int], numpy.ndarray], n: int):
    w = f(n)
    for i in range(n):
        print("{:6s} {:24s} {:6d} {:6d} {:53.40f}".format(source, name, n, i+1, w[i]))


def main():

    source = "numpy"
    filename = "numpy_windows.txt"
    with open(filename, "w") as fo, contextlib.redirect_stdout(fo):
        print("# Python: {}".format(sys.version))
        print("# Numpy: {}".format(numpy.__version__))
        for n in range(1, 100 + 1):
            dump_window_function(source, "rectwin"    , lambda n: numpy.ones(n)       , n)
            dump_window_function(source, "bartlett"   , lambda n: numpy.bartlett(n)   , n)
            dump_window_function(source, "blackman"   , lambda n: numpy.blackman(n)   , n)
            dump_window_function(source, "hamming"    , lambda n: numpy.hamming(n)    , n)
            dump_window_function(source, "hanning"    , lambda n: numpy.hanning(n)    , n)
            dump_window_function(source, "kaiser_0p5" , lambda n: numpy.kaiser(n, 0.5), n)
            dump_window_function(source, "kaiser_0p8" , lambda n: numpy.kaiser(n, 0.8), n)

    source = "scipy"
    filename = "scipy_windows.txt"
    with open(filename, "w") as fo, contextlib.redirect_stdout(fo):
        print("# Python: {}".format(sys.version))
        print("# Scipy: {}".format(scipy.__version__))
        for n in range(1, 100 + 1):
            dump_window_function(source, "barthann"                 , lambda n: scipy.signal.windows.barthann(n       ), n)
            dump_window_function(source, "barthann_periodic"        , lambda n: scipy.signal.windows.barthann(n, False), n)
            dump_window_function(source, "barthann_symmetric"       , lambda n: scipy.signal.windows.barthann(n, True ), n)

            dump_window_function(source, "bartlett"                 , lambda n: scipy.signal.windows.bartlett(n       ), n)
            dump_window_function(source, "bartlett_periodic"        , lambda n: scipy.signal.windows.bartlett(n, False), n)
            dump_window_function(source, "bartlett_symmetric"       , lambda n: scipy.signal.windows.bartlett(n, True ), n)

            dump_window_function(source, "blackman"                 , lambda n: scipy.signal.windows.blackman(n       ), n)
            dump_window_function(source, "blackman_periodic"        , lambda n: scipy.signal.windows.blackman(n, False), n)
            dump_window_function(source, "blackman_symmetric"       , lambda n: scipy.signal.windows.blackman(n, True ), n)

            dump_window_function(source, "blackmanharris"           , lambda n: scipy.signal.windows.blackmanharris(n       ), n)
            dump_window_function(source, "blackmanharris_periodic"  , lambda n: scipy.signal.windows.blackmanharris(n, False), n)
            dump_window_function(source, "blackmanharris_symmetric" , lambda n: scipy.signal.windows.blackmanharris(n, True ), n)

            dump_window_function(source, "bohman"                   , lambda n: scipy.signal.windows.bohman(n       ), n)
            dump_window_function(source, "bohman_periodic"          , lambda n: scipy.signal.windows.bohman(n, False), n)
            dump_window_function(source, "bohman_symmetric"         , lambda n: scipy.signal.windows.bohman(n, True ), n)

            dump_window_function(source, "boxcar"                   , lambda n: scipy.signal.windows.boxcar(n       ), n)
            dump_window_function(source, "boxcar_periodic"          , lambda n: scipy.signal.windows.boxcar(n, False), n)
            dump_window_function(source, "boxcar_symmetric"         , lambda n: scipy.signal.windows.boxcar(n, True ), n)

            dump_window_function(source, "chebwin_100p0"            , lambda n: scipy.signal.windows.chebwin(n, 100.0), n)
            dump_window_function(source, "chebwin_120p0"            , lambda n: scipy.signal.windows.chebwin(n, 120.0), n)

            dump_window_function(source, "chebwin_100p0_periodic"   , lambda n: scipy.signal.windows.chebwin(n, 100.0, False), n)
            dump_window_function(source, "chebwin_120p0_periodic"   , lambda n: scipy.signal.windows.chebwin(n, 120.0, False), n)

            dump_window_function(source, "chebwin_100p0_symmetric"  , lambda n: scipy.signal.windows.chebwin(n, 100.0, True), n)
            dump_window_function(source, "chebwin_120p0_symmetric"  , lambda n: scipy.signal.windows.chebwin(n, 120.0, True), n)

if __name__ == "__main__":
    main()
