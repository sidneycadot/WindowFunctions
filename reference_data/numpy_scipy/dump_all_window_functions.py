#! /usr/bin/env python3

import sys
import contextlib
import datetime
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
        print("# Python version {}".format(sys.version))
        print("# Numpy version {}".format(numpy.__version__))
        print("# {}".format(datetime.datetime.now().isoformat(' ')))
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
        print("# Python version {}".format(sys.version))
        print("# Numpy version {}".format(numpy.__version__))
        print("# Scipy version {}".format(scipy.__version__))
        print("# {}".format(datetime.datetime.now().isoformat(' ')))
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
            dump_window_function(source, "chebwin_100p0_periodic"   , lambda n: scipy.signal.windows.chebwin(n, 100.0, False), n)
            dump_window_function(source, "chebwin_100p0_symmetric"  , lambda n: scipy.signal.windows.chebwin(n, 100.0, True), n)

            dump_window_function(source, "chebwin_120p0"            , lambda n: scipy.signal.windows.chebwin(n, 120.0), n)
            dump_window_function(source, "chebwin_120p0_periodic"   , lambda n: scipy.signal.windows.chebwin(n, 120.0, False), n)
            dump_window_function(source, "chebwin_120p0_symmetric"  , lambda n: scipy.signal.windows.chebwin(n, 120.0, True), n)

            # TODO: dump values for scipy.signal.window.cosine() function.
            # TODO: dump values for scipy.signal.window.dpss() function.
            # TODO: dump values for scipy.signal.window.exponential() function.

            dump_window_function(source, "flattop"                  , lambda n: scipy.signal.windows.flattop(n       ), n)
            dump_window_function(source, "flattop_periodic"         , lambda n: scipy.signal.windows.flattop(n, False), n)
            dump_window_function(source, "flattop_symmetric"        , lambda n: scipy.signal.windows.flattop(n, True ), n)

            dump_window_function(source, "gausswin_2p5"             , lambda n: scipy.signal.windows.gaussian(n, 2.5       ), n)
            dump_window_function(source, "gausswin_2p5_periodic"    , lambda n: scipy.signal.windows.gaussian(n, 2.5, False), n)
            dump_window_function(source, "gausswin_2p5_symmetric"   , lambda n: scipy.signal.windows.gaussian(n, 2.5, True ), n)

            dump_window_function(source, "gausswin_3p2"             , lambda n: scipy.signal.windows.gaussian(n, 3.2       ), n)
            dump_window_function(source, "gausswin_3p2_periodic"    , lambda n: scipy.signal.windows.gaussian(n, 3.2, False), n)
            dump_window_function(source, "gausswin_3p2_symmetric"   , lambda n: scipy.signal.windows.gaussian(n, 3.2, True ), n)

            # TODO: dump values for scipy.signal.window.general_cosine() function.
            # TODO: dump values for scipy.signal.window.general_gaussian() function.
            # TODO: dump values for scipy.signal.window.general_hamming() function.

            dump_window_function(source, "hamming"               , lambda n: scipy.signal.windows.hamming(n       ), n)
            dump_window_function(source, "hamming_periodic"      , lambda n: scipy.signal.windows.hamming(n, False), n)
            dump_window_function(source, "hamming_symmetric"     , lambda n: scipy.signal.windows.hamming(n, True ), n)

            dump_window_function(source, "hann"                  , lambda n: scipy.signal.windows.hann(n       ), n)
            dump_window_function(source, "hann_periodic"         , lambda n: scipy.signal.windows.hann(n, False), n)
            dump_window_function(source, "hann_symmetric"        , lambda n: scipy.signal.windows.hann(n, True ), n)

            dump_window_function(source, "kaiser_0p5"            , lambda n: scipy.signal.windows.kaiser(n, 0.5       ), n)
            dump_window_function(source, "kaiser_0p5_periodic"   , lambda n: scipy.signal.windows.kaiser(n, 0.5, False), n)
            dump_window_function(source, "kaiser_0p5_symmetric"  , lambda n: scipy.signal.windows.kaiser(n, 0.5, True ), n)

            dump_window_function(source, "kaiser_0p8"            , lambda n: scipy.signal.windows.kaiser(n, 0.8       ), n)
            dump_window_function(source, "kaiser_0p8_periodic"   , lambda n: scipy.signal.windows.kaiser(n, 0.8, False), n)
            dump_window_function(source, "kaiser_0p8_symmetric"  , lambda n: scipy.signal.windows.kaiser(n, 0.8, True ), n)

            # TODO: dump values for scipy.signal.window.lanczos() function.

            dump_window_function(source, "nuttall"               , lambda n: scipy.signal.windows.nuttall(n       ), n)
            dump_window_function(source, "nuttall_periodic"      , lambda n: scipy.signal.windows.nuttall(n, False), n)
            dump_window_function(source, "nuttall_symmetric"     , lambda n: scipy.signal.windows.nuttall(n, True ), n)

            dump_window_function(source, "parzen"                , lambda n: scipy.signal.windows.parzen (n       ), n)
            dump_window_function(source, "parzen_periodic"       , lambda n: scipy.signal.windows.parzen (n, False), n)
            dump_window_function(source, "parzen_symmetric"      , lambda n: scipy.signal.windows.parzen (n, True ), n)

            #  TODO: dump values for scipy.signal.window.taylor() function.

            dump_window_function(source, "triang"                , lambda n: scipy.signal.windows.triang(n       ), n)
            dump_window_function(source, "triang_periodic"       , lambda n: scipy.signal.windows.triang(n, False), n)
            dump_window_function(source, "triang_symmetric"      , lambda n: scipy.signal.windows.triang(n, True ), n)

            dump_window_function(source, "tukey"                 , lambda n: scipy.signal.windows.tukey(n            ), n)

            dump_window_function(source, "tukey_0p0"             , lambda n: scipy.signal.windows.tukey(n, 0.0       ), n)
            dump_window_function(source, "tukey_0p0_periodic"    , lambda n: scipy.signal.windows.tukey(n, 0.0, False), n)
            dump_window_function(source, "tukey_0p0_symmetric"   , lambda n: scipy.signal.windows.tukey(n, 0.0, True ), n)

            dump_window_function(source, "tukey_0p2"             , lambda n: scipy.signal.windows.tukey(n, 0.2       ), n)
            dump_window_function(source, "tukey_0p2_periodic"    , lambda n: scipy.signal.windows.tukey(n, 0.2, False), n)
            dump_window_function(source, "tukey_0p2_symmetric"   , lambda n: scipy.signal.windows.tukey(n, 0.2, True ), n)

            dump_window_function(source, "tukey_0p5"             , lambda n: scipy.signal.windows.tukey(n, 0.5       ), n)
            dump_window_function(source, "tukey_0p5_periodic"    , lambda n: scipy.signal.windows.tukey(n, 0.5, False), n)
            dump_window_function(source, "tukey_0p5_symmetric"   , lambda n: scipy.signal.windows.tukey(n, 0.5, True ), n)

            dump_window_function(source, "tukey_0p8"             , lambda n: scipy.signal.windows.tukey(n, 0.8       ), n)
            dump_window_function(source, "tukey_0p8_periodic"    , lambda n: scipy.signal.windows.tukey(n, 0.8, False), n)
            dump_window_function(source, "tukey_0p8_symmetric"   , lambda n: scipy.signal.windows.tukey(n, 0.8, True ), n)

            dump_window_function(source, "tukey_1p0"             , lambda n: scipy.signal.windows.tukey(n, 1.0       ), n)
            dump_window_function(source, "tukey_1p0_periodic"    , lambda n: scipy.signal.windows.tukey(n, 1.0, False), n)
            dump_window_function(source, "tukey_1p0_symmetric"   , lambda n: scipy.signal.windows.tukey(n, 1.0, True ), n)


if __name__ == "__main__":
    main()
