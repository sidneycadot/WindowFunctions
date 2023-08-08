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

        # Check that this window is new.
        if not np.isnan(window[i]):
            raise RuntimeError("Duplicate reference value.")

        window[i] = value


def read_window_function_reference_data(filename: str) -> dict[str, WindowFunctionReferenceData]:
    """Read window-function reference data from a file."""

    reference_windows = {}

    for line in open(filename):

        (name, M, i, value) = line.split()
        M = int(M)
        i = int(i)
        value = float(value)

        if name in reference_windows:
            reference_window = reference_windows[name]
        else:
            reference_window = WindowFunctionReferenceData()
            reference_windows[name] = reference_window

        reference_window.set_value(M, i, value)

    # Done reading data.
    # Verify that all reference winfunc values are finite (not NaN).

    for reference_window in reference_windows.values():
        for w in reference_window.windows.values():
            assert np.isfinite(w).all()

    # Return the reference window functions read in.
    return reference_windows


def make_python_function_map_for_matlab() -> dict[str, Callable[[int], np.ndarray]]:
    """Make Python function map for Matlab window-functions.
    
    This lists Python equivalents for the window types found in the Matlab window-function reference data file.
    """

    winfunc_map = {}

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
    winfunc_map[ "hanning"                  ] = lambda M : wf.hanning(M)
    winfunc_map[ "hanning_periodic"         ] = lambda M : wf.hanning(M, False)
    winfunc_map[ "hanning_symmetric"        ] = lambda M : wf.hanning(M, True)
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


def make_python_function_map_for_octave()-> dict[str, Callable[[int], np.ndarray]]:
    """Make Python function map for Octave window-functions.
    
    This lists Python equivalents for the window types found in the Octave window-function reference data file.

    Differences with the Matlab functions:

    * The flattop window is defined differently (a cosine window with slightly different doefficients).
    * The Hanning window behaves identically to the Hann  window in Octave.
    * The Nutall window is defined differently (a cosine window with slightly different doefficients).
    * The Taylor window is not implemented in Octave.

    Note: In 2017, the Gauss window behaved differently in Octave compared to the Matlan version.
          This has since been fixed.
    """

    winfunc_map = {}

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
    winfunc_map[ "gausswin"                 ] = lambda M : wf.gausswin(M, 2.5)
    winfunc_map[ "gausswin_2p5"             ] = lambda M : wf.gausswin(M, 2.5)
    winfunc_map[ "gausswin_3p2"             ] = lambda M : wf.gausswin(M, 3.2)
    winfunc_map[ "hamming"                  ] = lambda M : wf.hamming(M)
    winfunc_map[ "hamming_periodic"         ] = lambda M : wf.hamming(M, False)
    winfunc_map[ "hamming_symmetric"        ] = lambda M : wf.hamming(M, True)
    winfunc_map[ "hanning"                  ] = lambda M : wf.hanning_octave(M)
    winfunc_map[ "hanning_periodic"         ] = lambda M : wf.hanning_octave(M, False)
    winfunc_map[ "hanning_symmetric"        ] = lambda M : wf.hanning_octave(M, True)
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


def check_reference_wf(reference_windows: dict[str, WindowFunctionReferenceData],
                       python_function_map: dict[str, Callable[[int], np.ndarray]]) -> None:
    """verify window-function reference data against corresponding Python implementations."""

    for (reference_window_name, reference_window) in reference_windows.items():

        if reference_window_name not in python_function_map:
            print("Not found in python_function_map ...... : {:30}".format(reference_window_name))
            continue

        python_function = python_function_map[reference_window_name]
        ok = True
        worsterr = 0.0
        for (M, reference_values) in reference_window.windows.items():
            python_values = python_function(M)
            err = np.abs(reference_values - python_values)
            maxerr = max(err) # Maximum error for this window size
            if maxerr > 1e-13:
                print("Bad window value: {} M = {} maxerr = {}".format(reference_window_name, M, maxerr))
                print("reference_values:", reference_values)
                print("python_values:", python_values)
                ok = False
            worsterr = max(maxerr, worsterr) # Maximum error across all window sizes.
        if ok:
            print("Perfect correspondence ........ : {:30} (worsterr = {:10.10g})".format(reference_window_name, worsterr))
    print()


def main() -> None:
    """Verify Matlab and Octave window-function reference data against Python implementations."""

    print()

    print("*** CHECK MATLAB REFERENCE VALUES ***")
    print()
    reference_windows = read_window_function_reference_data("reference/matlab_windows.txt")
    python_function_map = make_python_function_map_for_matlab()
    check_reference_wf(reference_windows, python_function_map)

    print("*** CHECK OCTAVE REFERENCE VALUES ***")
    print()
    reference_windows = read_window_function_reference_data("reference/octave_windows.txt")
    python_function_map = make_python_function_map_for_octave()
    check_reference_wf(reference_windows, python_function_map)


if __name__ == "__main__":
    main()
