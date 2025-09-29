#! /usr/bin/env -S python3 -B

"""Test Matlab, Octave, Python-NumPy, and Python-SciPy window-function reference data against Python implementations."""

import glob

import numpy as np
import window_functions as wf


class WindowFunctionReferenceData:
    """A reference window function contains reference data for different window sizes."""

    def __init__(self):
        self.windows: dict[int, np.ndarray] = {}

    def set_value(self, M: int , i: int, value: float) -> None:
        """Set single value inside an M-size window."""

        if M in self.windows:
            window = self.windows[M]
        else:
            window = np.full(M, np.nan) # initialize window with NaN values.
            self.windows[M] = window

        i -= 1 # Matlab/octave use 1-based indexing, but we use 0-based indexing.

        if not (0 <= i < M):
            raise ValueError("Bad reference value.")

        # Check that this data point is new.
        if not np.isnan(window[i]):
            raise ValueError("Duplicate reference value.")

        window[i] = value


def read_window_function_reference_data(filename: str) -> dict[str, WindowFunctionReferenceData]:
    """Read window-function reference data from a file."""

    reference_windows = {}

    for line in open(filename):

        if line.startswith("#"):
            # Ignore comment lines.
            continue

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

        try:
            reference_window.set_value(M, i, value)
        except ValueError:
            print(f"Error processing line: {line!r}")
            raise

    # Done reading data.
    # Verify that all reference window values are finite (not NaN).

    for reference_window in reference_windows.values():
        for w in reference_window.windows.values():
            assert np.isfinite(w).all()

    # Return the reference window functions read in.
    return reference_windows


_python_function_map = {

    # Python equivalents for the window types in the window-function reference data.

    ("matlab", "barthannwin")               : lambda M : wf.barthannwin(M),
    ("matlab", "bartlett")                  : lambda M : wf.bartlett(M),
    ("matlab", "blackman")                  : lambda M : wf.blackman(M),
    ("matlab", "blackman_periodic")         : lambda M : wf.blackman(M, False),
    ("matlab", "blackman_symmetric")        : lambda M : wf.blackman(M, True),
    ("matlab", "blackmanharris")            : lambda M : wf.blackmanharris(M),
    ("matlab", "blackmanharris_periodic")   : lambda M : wf.blackmanharris(M, False),
    ("matlab", "blackmanharris_symmetric")  : lambda M : wf.blackmanharris(M, True),
    ("matlab", "bohmanwin")                 : lambda M : wf.bohmanwin(M),
    ("matlab", "boxcar")                    : lambda M : wf.rectwin(M),
    ("matlab", "chebwin")                   : lambda M : wf.chebwin(M),
    ("matlab", "chebwin_100p0")             : lambda M : wf.chebwin(M, 100.0),
    ("matlab", "chebwin_120p0")             : lambda M : wf.chebwin(M, 120.0),
    ("matlab", "flattopwin")                : lambda M : wf.flattopwin(M),
    ("matlab", "flattopwin_periodic")       : lambda M : wf.flattopwin(M, False),
    ("matlab", "flattopwin_symmetric")      : lambda M : wf.flattopwin(M, True),
    ("matlab", "gausswin")                  : lambda M : wf.gausswin(M),
    ("matlab", "gausswin_2p5")              : lambda M : wf.gausswin(M, 2.5),
    ("matlab", "gausswin_3p2")              : lambda M : wf.gausswin(M, 3.2),
    ("matlab", "hamming")                   : lambda M : wf.hamming(M),
    ("matlab", "hamming_periodic")          : lambda M : wf.hamming(M, False),
    ("matlab", "hamming_symmetric")         : lambda M : wf.hamming(M, True),
    ("matlab", "hann")                      : lambda M : wf.hann(M),
    ("matlab", "hann_periodic")             : lambda M : wf.hann(M, False),
    ("matlab", "hann_symmetric")            : lambda M : wf.hann(M, True),
    ("matlab", "hanning")                   : lambda M : wf.hanning(M),
    ("matlab", "hanning_periodic")          : lambda M : wf.hanning(M, False),
    ("matlab", "hanning_symmetric")         : lambda M : wf.hanning(M, True),
    ("matlab", "kaiser")                    : lambda M : wf.kaiser(M),
    ("matlab", "kaiser_0p5")                : lambda M : wf.kaiser(M, 0.5),
    ("matlab", "kaiser_0p8")                : lambda M : wf.kaiser(M, 0.8),
    ("matlab", "nuttallwin")                : lambda M : wf.nuttallwin(M),
    ("matlab", "nuttallwin_periodic")       : lambda M : wf.nuttallwin(M, False),
    ("matlab", "nuttallwin_symmetric")      : lambda M : wf.nuttallwin(M, True),
    ("matlab", "parzenwin")                 : lambda M : wf.parzenwin(M),
    ("matlab", "rectwin")                   : lambda M : wf.rectwin(M),
    ("matlab", "taylorwin")                 : lambda M : wf.taylorwin(M),
    ("matlab", "taylorwin_3")               : lambda M : wf.taylorwin(M, 3),
    ("matlab", "taylorwin_4")               : lambda M : wf.taylorwin(M, 4),
    ("matlab", "taylorwin_5")               : lambda M : wf.taylorwin(M, 5),
    ("matlab", "taylorwin_3_m20")           : lambda M : wf.taylorwin(M, 3, -20),
    ("matlab", "taylorwin_4_m20")           : lambda M : wf.taylorwin(M, 4, -20),
    ("matlab", "taylorwin_5_m20")           : lambda M : wf.taylorwin(M, 5, -20),
    ("matlab", "taylorwin_3_m30")           : lambda M : wf.taylorwin(M, 3, -30),
    ("matlab", "taylorwin_4_m30")           : lambda M : wf.taylorwin(M, 4, -30),
    ("matlab", "taylorwin_5_m30")           : lambda M : wf.taylorwin(M, 5, -30),
    ("matlab", "taylorwin_3_m40")           : lambda M : wf.taylorwin(M, 3, -40),
    ("matlab", "taylorwin_4_m40")           : lambda M : wf.taylorwin(M, 4, -40),
    ("matlab", "taylorwin_5_m40")           : lambda M : wf.taylorwin(M, 5, -40),
    ("matlab", "triang")                    : lambda M : wf.triang(M),
    ("matlab", "tukeywin")                  : lambda M : wf.tukeywin(M),
    ("matlab", "tukeywin_0p0")              : lambda M : wf.tukeywin(M, 0.0),
    ("matlab", "tukeywin_0p2")              : lambda M : wf.tukeywin(M, 0.2),
    ("matlab", "tukeywin_0p5")              : lambda M : wf.tukeywin(M, 0.5),
    ("matlab", "tukeywin_0p8")              : lambda M : wf.tukeywin(M, 0.8),
    ("matlab", "tukeywin_1p0")              : lambda M : wf.tukeywin(M, 1.0),

    # Python equivalents for the window types in the Octave window-function reference data.

    # Differences with the Matlab functions:

    # * The flattop window is defined differently (a cosine window with slightly different doefficients).
    # * The Hanning window behaves identically to the Hann window in Octave.
    # * The Nutall window is defined differently (a cosine window with slightly different doefficients).

    ("octave", "barthannwin")               : lambda M : wf.barthannwin(M),
    ("octave", "bartlett")                  : lambda M : wf.bartlett(M),
    ("octave", "blackman")                  : lambda M : wf.blackman(M),
    ("octave", "blackman_periodic")         : lambda M : wf.blackman(M, False),
    ("octave", "blackman_symmetric")        : lambda M : wf.blackman(M, True),
    ("octave", "blackmanharris")            : lambda M : wf.blackmanharris(M),
    ("octave", "blackmanharris_periodic")   : lambda M : wf.blackmanharris(M, False),
    ("octave", "blackmanharris_symmetric")  : lambda M : wf.blackmanharris(M, True),
    ("octave", "blackmannuttall")           : lambda M : wf.nuttallwin(M),              # Replicates Matlabs's nutall() function.
    ("octave", "blackmannuttall_periodic")  : lambda M : wf.nuttallwin(M, False),
    ("octave", "blackmannuttall_symmetric") : lambda M : wf.nuttallwin(M, True),
    ("octave", "bohmanwin")                 : lambda M : wf.bohmanwin(M),
    ("octave", "boxcar")                    : lambda M : wf.rectwin(M),
    ("octave", "chebwin")                   : lambda M : wf.chebwin(M),
    ("octave", "chebwin_100p0")             : lambda M : wf.chebwin(M, 100.0),
    ("octave", "chebwin_120p0")             : lambda M : wf.chebwin(M, 120.0),
    ("octave", "expwin_2p5")                : lambda M : wf.expwin(M, 2.5),
    ("octave", "expwin_3p2")                : lambda M : wf.expwin(M, 3.2),
    ("octave", "expwin_minus_10")           : lambda M : wf.expwin(M, -10.0),
    ("octave", "expwin_minus_20")           : lambda M : wf.expwin(M, -20.0),
    ("octave", "expwin_2p5_canonical")      : lambda M : wf.expwin(M, 2.5, True),
    ("octave", "expwin_3p2_canonical")      : lambda M : wf.expwin(M, 3.2, True),
    ("octave", "expwin_minus_10_canonical") : lambda M : wf.expwin(M, -10.0, True),
    ("octave", "expwin_minus_20_canonical") : lambda M : wf.expwin(M, -20.0, True),
    ("octave", "flattopwin")                : lambda M : wf.flattopwin_octave(M),
    ("octave", "flattopwin_periodic")       : lambda M : wf.flattopwin_octave(M, False),
    ("octave", "flattopwin_symmetric")      : lambda M : wf.flattopwin_octave(M, True),
    ("octave", "gaussian")                  : lambda M : wf.gaussian(M),
    ("octave", "gaussian_1p0")              : lambda M : wf.gaussian(M, 1.0),
    ("octave", "gaussian_2p5")              : lambda M : wf.gaussian(M, 2.5),
    ("octave", "gausswin")                  : lambda M : wf.gausswin(M),
    ("octave", "gausswin_2p5")              : lambda M : wf.gausswin(M, 2.5),
    ("octave", "gausswin_3p2")              : lambda M : wf.gausswin(M, 3.2),
    ("octave", "hamming")                   : lambda M : wf.hamming(M),
    ("octave", "hamming_periodic")          : lambda M : wf.hamming(M, False),
    ("octave", "hamming_symmetric")         : lambda M : wf.hamming(M, True),
    ("octave", "hanning")                   : lambda M : wf.hanning_octave(M),
    ("octave", "hanning_periodic")          : lambda M : wf.hanning_octave(M, False),
    ("octave", "hanning_symmetric")         : lambda M : wf.hanning_octave(M, True),
    ("octave", "hann")                      : lambda M : wf.hann(M),
    ("octave", "hann_periodic")             : lambda M : wf.hann(M, False),
    ("octave", "hann_symmetric")            : lambda M : wf.hann(M, True),
    ("octave", "kaiser")                    : lambda M : wf.kaiser(M),
    ("octave", "kaiser_0p5")                : lambda M : wf.kaiser(M, 0.5),
    ("octave", "kaiser_0p8")                : lambda M : wf.kaiser(M, 0.8),
    ("octave", "nuttallwin")                : lambda M : wf.nuttallwin_octave(M),
    ("octave", "nuttallwin_periodic")       : lambda M : wf.nuttallwin_octave(M, False),
    ("octave", "nuttallwin_symmetric")      : lambda M : wf.nuttallwin_octave(M, True),
    ("octave", "parzenwin")                 : lambda M : wf.parzenwin(M),
    ("octave", "poisswin_2p5")              : lambda M : wf.poisswin(M, 2.5),
    ("octave", "poisswin_3p2")              : lambda M : wf.poisswin(M, 3.2),
    ("octave", "rectwin")                   : lambda M : wf.rectwin(M),
    ("octave", "triang")                    : lambda M : wf.triang(M),
    ("octave", "taylorwin")                 : lambda M : wf.taylorwin(M),
    ("octave", "taylorwin_3")               : lambda M : wf.taylorwin(M, 3),
    ("octave", "taylorwin_4")               : lambda M : wf.taylorwin(M, 4),
    ("octave", "taylorwin_5")               : lambda M : wf.taylorwin(M, 5),
    ("octave", "taylorwin_3_m20")           : lambda M : wf.taylorwin(M, 3, -20),
    ("octave", "taylorwin_4_m20")           : lambda M : wf.taylorwin(M, 4, -20),
    ("octave", "taylorwin_5_m20")           : lambda M : wf.taylorwin(M, 5, -20),
    ("octave", "taylorwin_3_m30")           : lambda M : wf.taylorwin(M, 3, -30),
    ("octave", "taylorwin_4_m30")           : lambda M : wf.taylorwin(M, 4, -30),
    ("octave", "taylorwin_5_m30")           : lambda M : wf.taylorwin(M, 5, -30),
    ("octave", "taylorwin_3_m40")           : lambda M : wf.taylorwin(M, 3, -40),
    ("octave", "taylorwin_4_m40")           : lambda M : wf.taylorwin(M, 4, -40),
    ("octave", "taylorwin_5_m40")           : lambda M : wf.taylorwin(M, 5, -40),
    ("octave", "tukeywin")                  : lambda M : wf.tukeywin(M),
    ("octave", "tukeywin_0p0")              : lambda M : wf.tukeywin(M, 0.0),
    ("octave", "tukeywin_0p2")              : lambda M : wf.tukeywin(M, 0.2),
    ("octave", "tukeywin_0p5")              : lambda M : wf.tukeywin(M, 0.5),
    ("octave", "tukeywin_0p8")              : lambda M : wf.tukeywin(M, 0.8),
    ("octave", "tukeywin_1p0")              : lambda M : wf.tukeywin(M, 1.0),

    # Numpy windows.

    ("numpy", "ones")                       : lambda M : wf.rectwin(M),
    ("numpy", "bartlett")                   : lambda M : wf.bartlett(M),
    ("numpy", "blackman")                   : lambda M : wf.blackman(M),
    ("numpy", "hamming")                    : lambda M : wf.hamming(M),
    ("numpy", "hanning")                    : lambda M : wf.hanning_octave(M),
    ("numpy", "kaiser_0p5")                 : lambda M : wf.kaiser(M, 0.5),
    ("numpy", "kaiser_0p8")                 : lambda M : wf.kaiser(M, 0.8),

    # Scipy windows.

    ("scipy", "barthann")                   : lambda M : wf.barthannwin(M, True ),
    ("scipy", "barthann_periodic")          : lambda M : wf.barthannwin(M, False),
    ("scipy", "barthann_symmetric")         : lambda M : wf.barthannwin(M, True ),

    ("scipy", "bartlett")                   : lambda M : wf.bartlett(M, True ),
    ("scipy", "bartlett_periodic")          : lambda M : wf.bartlett(M, False),
    ("scipy", "bartlett_symmetric")         : lambda M : wf.bartlett(M, True ),

    ("scipy", "blackman")                   : lambda M : wf.blackman(M, True ),
    ("scipy", "blackman_periodic")          : lambda M : wf.blackman(M, False),
    ("scipy", "blackman_symmetric")         : lambda M : wf.blackman(M, True ),

    ("scipy", "blackmanharris")             : lambda M : wf.blackmanharris(M, True ),
    ("scipy", "blackmanharris_periodic")    : lambda M : wf.blackmanharris(M, False),
    ("scipy", "blackmanharris_symmetric")   : lambda M : wf.blackmanharris(M, True ),

    ("scipy", "bohman")                     : lambda M : wf.bohmanwin(M, True ),
    ("scipy", "bohman_periodic")            : lambda M : wf.bohmanwin(M, False),
    ("scipy", "bohman_symmetric")           : lambda M : wf.bohmanwin(M, True ),

    ("scipy", "boxcar")                     : lambda M : wf.rectwin(M),
    ("scipy", "boxcar_periodic")            : lambda M : wf.rectwin(M),
    ("scipy", "boxcar_symmetric")           : lambda M : wf.rectwin(M),

    ("scipy", "chebwin_100p0")              : lambda M : wf.chebwin(M, 100.0),
    ("scipy", "chebwin_100p0_periodic")     : lambda M : wf.chebwin(M, 100.0, False),
    ("scipy", "chebwin_100p0_symmetric")    : lambda M : wf.chebwin(M, 100.0, True),

    ("scipy", "chebwin_120p0")              : lambda M : wf.chebwin(M, 120.0),
    ("scipy", "chebwin_120p0_periodic")     : lambda M : wf.chebwin(M, 120.0, False),
    ("scipy", "chebwin_120p0_symmetric")    : lambda M : wf.chebwin(M, 120.0, True),

    ("scipy", "cosine")                     : lambda M : wf.cosine_scipy(M),
    ("scipy", "cosine_periodic")            : lambda M : wf.cosine_scipy(M, False),
    ("scipy", "cosine_symmetric")           : lambda M : wf.cosine_scipy(M, True),

    ("scipy", "lanczos")                     : lambda M : wf.lanczos(M),
    ("scipy", "lanczos_periodic")            : lambda M : wf.lanczos(M, False),
    ("scipy", "lanczos_symmetric")           : lambda M : wf.lanczos(M, True),

    ("scipy", "flattop")                    : lambda M : wf.flattopwin(M, True),
    ("scipy", "flattop_periodic")           : lambda M : wf.flattopwin(M, False),
    ("scipy", "flattop_symmetric")          : lambda M : wf.flattopwin(M, True),

    ("scipy", "gaussian_2p5")               : lambda M : wf.gaussian_scipy(M, 2.5),
    ("scipy", "gaussian_2p5_periodic")      : lambda M : wf.gaussian_scipy(M, 2.5, False),
    ("scipy", "gaussian_2p5_symmetric")     : lambda M : wf.gaussian_scipy(M, 2.5, True),
    ("scipy", "gaussian_3p2")               : lambda M : wf.gaussian_scipy(M, 3.2),
    ("scipy", "gaussian_3p2_periodic")      : lambda M : wf.gaussian_scipy(M, 3.2, False),
    ("scipy", "gaussian_3p2_symmetric")     : lambda M : wf.gaussian_scipy(M, 3.2, True),

    ("scipy", "hamming")                    : lambda M : wf.hamming(M, True),
    ("scipy", "hamming_periodic")           : lambda M : wf.hamming(M, False),
    ("scipy", "hamming_symmetric")          : lambda M : wf.hamming(M, True),

    ("scipy", "general_hamming_0p3")           : lambda M : wf.general_hamming_scipy(M, 0.3, True),
    ("scipy", "general_hamming_0p3_periodic")  : lambda M : wf.general_hamming_scipy(M, 0.3, False),
    ("scipy", "general_hamming_0p3_symmetric") : lambda M : wf.general_hamming_scipy(M, 0.3, True),
    ("scipy", "general_hamming_0p8")           : lambda M : wf.general_hamming_scipy(M, 0.8, True),
    ("scipy", "general_hamming_0p8_periodic")  : lambda M : wf.general_hamming_scipy(M, 0.8, False),
    ("scipy", "general_hamming_0p8_symmetric") : lambda M : wf.general_hamming_scipy(M, 0.8, True),

    ("scipy", "hann")                       : lambda M : wf.hann(M, True),
    ("scipy", "hann_periodic")              : lambda M : wf.hann(M, False),
    ("scipy", "hann_symmetric")             : lambda M : wf.hann(M, True),

    ("scipy", "kaiser_0p5")                 : lambda M : wf.kaiser(M, 0.5, True),
    ("scipy", "kaiser_0p5_periodic")        : lambda M : wf.kaiser(M, 0.5, False),
    ("scipy", "kaiser_0p5_symmetric")       : lambda M : wf.kaiser(M, 0.5, True),
    ("scipy", "kaiser_0p8")                 : lambda M : wf.kaiser(M, 0.8, True),
    ("scipy", "kaiser_0p8_periodic")        : lambda M : wf.kaiser(M, 0.8, False),
    ("scipy", "kaiser_0p8_symmetric")       : lambda M : wf.kaiser(M, 0.8, True),

    ("scipy", "nuttall")                    : lambda M : wf.nuttallwin(M, True ),
    ("scipy", "nuttall_periodic")           : lambda M : wf.nuttallwin(M, False),
    ("scipy", "nuttall_symmetric")          : lambda M : wf.nuttallwin(M, True ),

    ("scipy", "parzen")                     : lambda M : wf.parzenwin(M, True ),
    ("scipy", "parzen_periodic")            : lambda M : wf.parzenwin(M, False),
    ("scipy", "parzen_symmetric")           : lambda M : wf.parzenwin(M, True ),

    ("scipy", "triang")                     : lambda M : wf.triang(M, True ),
    ("scipy", "triang_periodic")            : lambda M : wf.triang(M, False),
    ("scipy", "triang_symmetric")           : lambda M : wf.triang(M, True ),

    ("scipy", "tukey")                      : lambda M : wf.tukeywin(M),
    ("scipy", "tukey_0p0")                  : lambda M : wf.tukeywin(M, 0.0),
    ("scipy", "tukey_0p0_periodic")         : lambda M : wf.tukeywin(M, 0.0, False),
    ("scipy", "tukey_0p0_symmetric")        : lambda M : wf.tukeywin(M, 0.0, True),
    ("scipy", "tukey_0p2")                  : lambda M : wf.tukeywin(M, 0.2),
    ("scipy", "tukey_0p2_periodic")         : lambda M : wf.tukeywin(M, 0.2, False),
    ("scipy", "tukey_0p2_symmetric")        : lambda M : wf.tukeywin(M, 0.2, True),
    ("scipy", "tukey_0p5")                  : lambda M : wf.tukeywin(M, 0.5),
    ("scipy", "tukey_0p5_periodic")         : lambda M : wf.tukeywin(M, 0.5, False),
    ("scipy", "tukey_0p5_symmetric")        : lambda M : wf.tukeywin(M, 0.5, True),
    ("scipy", "tukey_0p8")                  : lambda M : wf.tukeywin(M, 0.8),
    ("scipy", "tukey_0p8_periodic")         : lambda M : wf.tukeywin(M, 0.8, False),
    ("scipy", "tukey_0p8_symmetric")        : lambda M : wf.tukeywin(M, 0.8, True),
    ("scipy", "tukey_1p0")                  : lambda M : wf.tukeywin(M, 1.0),
    ("scipy", "tukey_1p0_periodic")         : lambda M : wf.tukeywin(M, 1.0, False),
    ("scipy", "tukey_1p0_symmetric")        : lambda M : wf.tukeywin(M, 1.0, True)
}


def check_reference_waveform_data(reference_windows: dict[[str, WindowFunctionReferenceData], np.ndarray]) -> None:
    """Verify window-function reference data against corresponding Python implementations."""

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
            print("Perfect correspondence : {:40} (worsterr = {:8.3g})".format("/".join(reference_window_key), worsterr))


def main() -> None:
    """Verify Window-function reference data against the Python implementations."""

    for filename in ["../reference_data/numpy_scipy/scipy_windows.txt"]:
    #for filename in glob.glob("reference_data/*_windows.txt"):
        print("Checking {} ...".format(filename))
        print()
        check_reference_waveform_data(read_window_function_reference_data(filename))
        print()


if __name__ == "__main__":
    main()
