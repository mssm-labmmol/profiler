#!/usr/bin/env python3

import numpy as np
from sys import exit

tol = 1e-03
ref = np.loadtxt('reference.dat')
comp = np.loadtxt('phe-ala-gen.dat')
err = np.abs((comp - ref))
for e in err:
    if e > tol:
        exit(1)

