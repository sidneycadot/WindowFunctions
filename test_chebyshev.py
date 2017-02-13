#! /usr/bin/env python3

import numpy as np
import window_functions as wf

print()
print("Replicating example values from Lyons, Appendix I:")
print()

w = wf.chebwin(9, 60)

for (i, x) in enumerate(w):
    print("{:3} --> {:.4f}".format(i, x))

print()

ref = np.array([0.0519, 0.2271, 0.5379, 0.8605, 1.0000, 0.8605, 0.5379, 0.2271, 0.0519])

print("reference value check:", (np.abs(ref-w) < 5e-5).all())
print()
