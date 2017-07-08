#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp2d

n = 1/np.array([10, 20])
N = len(n)

xb = .1
yb = .1 

fmms = lambda x, y: np.sin(np.pi*x/xb) * np.sin(np.pi*y/yb)

f_err = lambda a, b: np.fabs(a - b)/b

err = np.zeros(N)

pref = 'run'

for i in range(N):

	X = np.loadtxt(pref+str(i)+'/X0')
	Y = np.loadtxt(pref+str(i)+'/Y0')
	df = np.loadtxt(pref+str(i)+'/out0')

	func = interp2d(X, Y, df)

	err[i] = f_err(func(xb/2, yb/2), fmms(xb/2, yb/2))

fit = np.polyfit(np.log(n), np.log(err), 1)
print(fit[0])

print(err)

plt.loglog(n, err, '-o')
plt.show()


