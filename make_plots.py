#! /usr/bin/env python3

import window_functions as wf
import matplotlib.pyplot as plt
import numpy as np

def plot_cosine_windows(n):

    plt.title("cosine windows")
    plt.plot(wf.rectwin(n), label = "rectwin")
    plt.plot(wf.hamming(n), label = "hamming")
    plt.plot(wf.hann(n), label = "hann")
    plt.plot(wf.blackman(n), label = "blackman")
    plt.plot(wf.nuttallwin(n), label = "nutallwin")
    plt.plot(wf.blackmanharris(n), label = "blackmanharris")
    plt.plot(wf.nuttallwin_octave(n), label = "nutallwin (octave)")
    plt.plot(wf.flattopwin(n), label = "flattopwin")
    plt.plot(wf.flattopwin_octave(n), label = "flattopwin (octave)")
    plt.ylim(None, 1.2)
    plt.legend()
    plt.show()

def plot_chebyshev_windows(n):

    plt.title("chebyshev windows")
    for r in [10, 20, 50, 100, 200, 500, 1000]:
        plt.plot(wf.chebwin(n, r), label = "r = {}".format(r))
    plt.legend()

def plot_kaiser_windows(n):

    plt.title("kaiser windows")
    for beta in [1, 2, 5, 10, 20, 50]:
        plt.plot(wf.kaiser(n, beta), label = "beta = {:.2f}".format(beta))
    plt.legend()

def plot_tukey_windows(n):

    plt.title("tukey windows")
    for r in [0.25, 0.50, 0.75, 1.0]:
        plt.plot(wf.tukeywin(n, r), label = "r = {:.2f}".format(r))
    plt.legend()

def plot_gaussian_windows(n):

    plt.title("gaussian windows")
    plt.plot(wf.parzenwin(n), lw = 5, label = "parzenwin")
    for alpha in [0.5, 1.0, 2.0, 3.3, 5.0, 10.0]:
        plt.plot(wf.gausswin(n, alpha), label = "alpha = {:.2f}".format(alpha))
    plt.legend()

n = 1000
plt.subplot(321)
plot_chebyshev_windows(n)
plt.subplot(322)
plot_kaiser_windows(n)
plt.subplot(323)
plot_tukey_windows(n)
plt.subplot(324)
plot_gaussian_windows(n)
plt.show()

