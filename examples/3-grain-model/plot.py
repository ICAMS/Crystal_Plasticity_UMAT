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
fnames = ['Job-1_sig_eps.csv']


# initialize data arrays
Nlc = len(fnames)
data = dict()

# read data
for nc, fname in enumerate(fnames):
    data[nc] = np.loadtxt(open(fname, "rb"), delimiter=",", skiprows=1)
#Ndat = len(data[0])
#time = np.linspace(1, Ndat+1, Ndat, endpoint=True)

E = (data[0][4, 13] - data[0][0, 13]) / (
    data[0][4, 1] - data[0][4, 7] - data[0][0, 1] + data[0][0, 7])
print(f"Estimated Young's modulus E = {E} GPa")
E = 171250.
print(f"True Young's modulus E = {E} GPa")

# plot stress-strain
col = ['r', 'b', 'orange', 'g']
for i, val in enumerate(data.values()):
    seq = fe.sig_eq_j2(val[:,13:19])
    peeq = fe.eps_eq(val[:,7:13])
    eteq = fe.eps_eq(val[:,1:7])
    plt.plot(peeq, seq, color=col[2*i], linewidth=1,
             label='Epl')
    plt.plot(eteq, seq, color=col[2*i+1], linewidth=1,
             label='Etot')
    plt.plot(peeq+seq/E, seq, ':k', linewidth=1,
             label='Epl + Eel')
    plt.ylabel('stress (MPa)')
    plt.xlabel('strain (.)')
    plt.legend()
    plt.show()

"""
# plot time curves
for i, val in enumerate(data.values()):
    plt.plot(val[:, 0], 100*seq/E, '-b', linewidth=1,
             label='SIG')
    plt.plot(val[:, 0], peeq, '-r', linewidth=1,
             label='Epl')
    plt.plot(val[:, 0], eteq, '-k', linewidth=1,
             label='Etot')
    plt.plot(val[:, 0], 100*(eteq-peeq), '-g', linewidth=1,
             label='100 Eel')
plt.ylabel('stress/E, strain')
plt.xlabel('time (a.u.)')
plt.legend()
plt.show()
"""


