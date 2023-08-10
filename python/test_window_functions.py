#! /usr/bin/env python3

"""Test Matlab and Octave window-function reference data against Python implementations."""

from typing import Callable

import numpy as np
import window_functions as wf


class WindowFunctionReferenceData:
    """A reference window function contains reference data for different window sizes."""

    def __init__(self):
        self.windows = {}

    def set_value(self, M: int , i: int, value: float) -> None:
        """Set single value inside an M-size window."""

        if M in self.windows:
            window = self.windows[M]
        else:
            window = np.full(M, np.nan) # initialize window with NaN values.
            self.windows[M] = window

        i -= 1 # Matlab/octave use 1-based indexing, but we use 0-based indexing.

        if not (0 <= i < M):
            raise RuntimeError("Bad reference value.")

        # Check that this data point is new.
        if not np.isnan(window[i]):
            raise RuntimeError("Duplicate reference value.")

        window[i] = value


def read_window_function_reference_data(filename: str) -> dict[str, WindowFunctionReferenceData]:
    """Read window-function reference data from a file."""

    reference_windows = {}

    for line in open(filename):

        (source, name, M, i, value) = line.split()
        M = int(M)
        i = int(i)
        value = float(value)

        key = (source, name)

        if key in reference_windows:
            reference_window = reference_windows[key]
        else:
            reference_window = WindowFunctionReferenceData()
            reference_windows[key] = reference_window

        reference_window.set_value(M, i, value)

    # Done reading data.
    # Verify that all reference winfunc values are finite (not NaN).

    for reference_window in reference_windows.values():
        for w in reference_window.windows.values():
            assert np.isfinite(w).all()

    # Return the reference window functions read in.
    return reference_windows


_python_function_map = {

    # Python equivalents for the window types in the Matlab window-function reference data.

    ("matlab", "barthannwin")              : lambda M : wf.barthannwin(M),
    ("matlab", "bartlett")                 : lambda M : wf.bartlett(M),
    ("matlab", "blackman")                 : lambda M : wf.blackman(M),
    ("matlab", "blackman_periodic")        : lambda M : wf.blackman(M, False),
    ("matlab", "blackman_symmetric")       : lambda M : wf.blackman(M, True),
    ("matlab", "blackmanharris")           : lambda M : wf.blackmanharris(M),
    ("matlab", "blackmanharris_periodic")  : lambda M : wf.blackmanharris(M, False),
    ("matlab", "blackmanharris_symmetric") : lambda M : wf.blackmanharris(M, True),
    ("matlab", "bohmanwin")                : lambda M : wf.bohmanwin(M),
    ("matlab", "chebwin")                  : lambda M : wf.chebwin(M),
    ("matlab", "chebwin_100p0")            : lambda M : wf.chebwin(M, 100.0),
    ("matlab", "chebwin_120p0")            : lambda M : wf.chebwin(M, 120.0),
    ("matlab", "flattopwin")               : lambda M : wf.flattopwin(M),
    ("matlab", "flattopwin_periodic")      : lambda M : wf.flattopwin(M, False),
    ("matlab", "flattopwin_symmetric")     : lambda M : wf.flattopwin(M, True),
    ("matlab", "gausswin")                 : lambda M : wf.gausswin(M),
    ("matlab", "gausswin_2p5")             : lambda M : wf.gausswin(M, 2.5),
    ("matlab", "gausswin_3p2")             : lambda M : wf.gausswin(M, 3.2),
    ("matlab", "hamming")                  : lambda M : wf.hamming(M),
    ("matlab", "hamming_periodic")         : lambda M : wf.hamming(M, False),
    ("matlab", "hamming_symmetric")        : lambda M : wf.hamming(M, True),
    ("matlab", "hann")                     : lambda M : wf.hann(M),
    ("matlab", "hann_periodic")            : lambda M : wf.hann(M, False),
    ("matlab", "hann_symmetric")           : lambda M : wf.hann(M, True),
    ("matlab", "hanning")                  : lambda M : wf.hanning(M),
    ("matlab", "hanning_periodic")         : lambda M : wf.hanning(M, False),
    ("matlab", "hanning_symmetric")        : lambda M : wf.hanning(M, True),
    ("matlab", "kaiser")                   : lambda M : wf.kaiser(M),
    ("matlab", "kaiser_0p5")               : lambda M : wf.kaiser(M, 0.5),
    ("matlab", "kaiser_0p8")               : lambda M : wf.kaiser(M, 0.8),
    ("matlab", "nuttallwin")               : lambda M : wf.nuttallwin(M),
    ("matlab", "nuttallwin_periodic")      : lambda M : wf.nuttallwin(M, False),
    ("matlab", "nuttallwin_symmetric")     : lambda M : wf.nuttallwin(M, True),
    ("matlab", "parzenwin")                : lambda M : wf.parzenwin(M),
    ("matlab", "rectwin")                  : lambda M : wf.rectwin(M),
    ("matlab", "taylorwin")                : lambda M : wf.taylorwin(M),
    ("matlab", "taylorwin_3")              : lambda M : wf.taylorwin(M, 3),
    ("matlab", "taylorwin_4")              : lambda M : wf.taylorwin(M, 4),
    ("matlab", "taylorwin_5")              : lambda M : wf.taylorwin(M, 5),
    ("matlab", "taylorwin_3_m20")          : lambda M : wf.taylorwin(M, 3, -20),
    ("matlab", "taylorwin_4_m20")          : lambda M : wf.taylorwin(M, 4, -20),
    ("matlab", "taylorwin_5_m20")          : lambda M : wf.taylorwin(M, 5, -20),
    ("matlab", "taylorwin_3_m30")          : lambda M : wf.taylorwin(M, 3, -30),
    ("matlab", "taylorwin_4_m30")          : lambda M : wf.taylorwin(M, 4, -30),
    ("matlab", "taylorwin_5_m30")          : lambda M : wf.taylorwin(M, 5, -30),
    ("matlab", "taylorwin_3_m40")          : lambda M : wf.taylorwin(M, 3, -40),
    ("matlab", "taylorwin_4_m40")          : lambda M : wf.taylorwin(M, 4, -40),
    ("matlab", "taylorwin_5_m40")          : lambda M : wf.taylorwin(M, 5, -40),
    ("matlab", "triang")                   : lambda M : wf.triang(M),
    ("matlab", "tukeywin")                 : lambda M : wf.tukeywin(M),
    ("matlab", "tukeywin_0p0")             : lambda M : wf.tukeywin(M, 0.0),
    ("matlab", "tukeywin_0p2")             : lambda M : wf.tukeywin(M, 0.2),
    ("matlab", "tukeywin_0p5")             : lambda M : wf.tukeywin(M, 0.5),
    ("matlab", "tukeywin_0p8")             : lambda M : wf.tukeywin(M, 0.8),
    ("matlab", "tukeywin_1p0")             : lambda M : wf.tukeywin(M, 1.0),

    # Python equivalents for the window types in the Octave window-function reference data.

    # Differences with the Matlab functions:

    # * The flattop window is defined differently (a cosine window with slightly different doefficients).
    # * The Hanning window behaves identically to the Hann  window in Octave.
    # * The Nutall window is defined differently (a cosine window with slightly different doefficients).
    # * The Taylor window is not implemented in Octave.

    # Note: In 2017, the Gauss window behaved differently in Octave compared to the Matlab version.
    #      This has since been fixed.

    ("octave", "barthannwin")              : lambda M : wf.barthannwin(M),
    ("octave", "bartlett")                 : lambda M : wf.bartlett(M),
    ("octave", "blackman")                 : lambda M : wf.blackman(M),
    ("octave", "blackman_periodic")        : lambda M : wf.blackman(M, False),
    ("octave", "blackman_symmetric")       : lambda M : wf.blackman(M, True),
    ("octave", "blackmanharris")           : lambda M : wf.blackmanharris(M),
    ("octave", "blackmanharris_periodic")  : lambda M : wf.blackmanharris(M, False),
    ("octave", "blackmanharris_symmetric") : lambda M : wf.blackmanharris(M, True),
    ("octave", "bohmanwin")                : lambda M : wf.bohmanwin(M),
    ("octave", "chebwin")                  : lambda M : wf.chebwin(M),
    ("octave", "chebwin_100p0")            : lambda M : wf.chebwin(M, 100.0),
    ("octave", "chebwin_120p0")            : lambda M : wf.chebwin(M, 120.0),
    ("octave", "flattopwin")               : lambda M : wf.flattopwin_octave(M),
    ("octave", "flattopwin_periodic")      : lambda M : wf.flattopwin_octave(M, False),
    ("octave", "flattopwin_symmetric")     : lambda M : wf.flattopwin_octave(M, True),
    ("octave", "gausswin")                 : lambda M : wf.gausswin(M, 2.5),
    ("octave", "gausswin_2p5")             : lambda M : wf.gausswin(M, 2.5),
    ("octave", "gausswin_3p2")             : lambda M : wf.gausswin(M, 3.2),
    ("octave", "hamming")                  : lambda M : wf.hamming(M),
    ("octave", "hamming_periodic")         : lambda M : wf.hamming(M, False),
    ("octave", "hamming_symmetric")        : lambda M : wf.hamming(M, True),
    ("octave", "hanning")                  : lambda M : wf.hanning_octave(M),
    ("octave", "hanning_periodic")         : lambda M : wf.hanning_octave(M, False),
    ("octave", "hanning_symmetric")        : lambda M : wf.hanning_octave(M, True),
    ("octave", "hann")                     : lambda M : wf.hann(M),
    ("octave", "hann_periodic")            : lambda M : wf.hann(M, False),
    ("octave", "hann_symmetric")           : lambda M : wf.hann(M, True),
    ("octave", "kaiser")                   : lambda M : wf.kaiser(M),
    ("octave", "kaiser_0p5")               : lambda M : wf.kaiser(M, 0.5),
    ("octave", "kaiser_0p8")               : lambda M : wf.kaiser(M, 0.8),
    ("octave", "nuttallwin")               : lambda M : wf.nuttallwin_octave(M),
    ("octave", "nuttallwin_periodic")      : lambda M : wf.nuttallwin_octave(M, False),
    ("octave", "nuttallwin_symmetric")     : lambda M : wf.nuttallwin_octave(M, True),
    ("octave", "parzenwin")                : lambda M : wf.parzenwin(M),
    ("octave", "rectwin")                  : lambda M : wf.rectwin(M),
    ("octave", "triang")                   : lambda M : wf.triang(M),
    ("octave", "tukeywin")                 : lambda M : wf.tukeywin(M),
    ("octave", "tukeywin_0p0")             : lambda M : wf.tukeywin(M, 0.0),
    ("octave", "tukeywin_0p2")             : lambda M : wf.tukeywin(M, 0.2),
    ("octave", "tukeywin_0p5")             : lambda M : wf.tukeywin(M, 0.5),
    ("octave", "tukeywin_0p8")             : lambda M : wf.tukeywin(M, 0.8),
    ("octave", "tukeywin_1p0")             : lambda M : wf.tukeywin(M, 1.0)
}


def check_reference_wf(reference_windows: dict[[str, WindowFunctionReferenceData], np.ndarray]) -> None:
    """verify window-function reference data against corresponding Python implementations."""

    for (reference_window_key, reference_window) in reference_windows.items():

        if reference_window_key not in _python_function_map:
            print("Not found in _python_function_map ...... : {}".format("/".join(reference_window_key)))
            continue

        python_function = _python_function_map[reference_window_key]
        ok = True
        worsterr = 0.0
        for (M, reference_values) in reference_window.windows.items():
            python_values = python_function(M)
            err = np.abs(reference_values - python_values)
            maxerr = max(err) # Maximum error for this window size
            if maxerr > 1e-13:
                print("Bad window value: {} M = {} maxerr = {}".format("/".join(reference_window_key), M, maxerr))
                print("reference_values:", reference_values)
                print("python_values:", python_values)
                ok = False
            worsterr = max(maxerr, worsterr) # Maximum error across all window sizes.
        if ok:
            print("Perfect correspondence ........ : {:30} (worsterr = {:10.10g})".format("/".join(reference_window_key), worsterr))
    print()


def main() -> None:
    """Verify Matlab and Octave window-function reference data against Python implementations."""

    print()

    #print("*** CHECK MATLAB REFERENCE VALUES ***")
    #print()
    #reference_windows = read_window_function_reference_data("reference_data/matlab_windows.txt")
    #python_function_map = make_python_function_map_for_matlab()
    #check_reference_wf(reference_windows, python_function_map)

    print("*** CHECK OCTAVE REFERENCE VALUES ***")
    print()
    reference_windows = read_window_function_reference_data("reference_data/octave_windows.txt")
    check_reference_wf(reference_windows)


if __name__ == "__main__":
    main()
