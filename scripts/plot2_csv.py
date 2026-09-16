#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot stress strain hysteresis loops

Created on Mon Aug  8 07:33:17 2022

@author: alexander
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import pylabfea as fe

# define list of filenames 
fnames = ['test1_glob.csv', 'test4_glob.csv']

# initialize data arrays
Nlc = len(fnames)
data = dict()

# read data
for nc, fname in enumerate(fnames):
    data[nc] = np.loadtxt(open(fname, "rb"), delimiter=",", skiprows=1)

# plot stress-strain
col = ['r', 'b', 'orange', 'g']
stl = ['-', '--']
for i, val in enumerate(data.values()):
    sig = val[:, 14]
    eps = val[:, 2]
    epl = val[:, 8]
    plt.plot(eps, sig, color=col[2*i], linewidth=1,
             linestyle=stl[i%2], label='Etot '+fnames[i][:5])
plt.ylabel('stress (MPa)')
plt.xlabel('strain (.)')
plt.legend()
plt.show()

for i, val in enumerate(data.values()):
    sig = val[:, 14]
    eps = val[:, 2]
    epl = val[:, 8]
    plt.plot(epl, sig, color=col[2*i+1], linewidth=1,
             linestyle=stl[i%2], label='Epl '+fnames[i][:5])
plt.ylabel('stress (MPa)')
plt.xlabel('strain (.)')
plt.legend()
plt.show()
