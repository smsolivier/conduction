#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import os 

N = 149

datadir = 'data/'

for i in range(N):

	print(i/N, end='\r')

	X = np.loadtxt('data/X'+str(i))
	Y = np.loadtxt('data/Y'+str(i))
	df = np.loadtxt('data/out'+str(i))

	plt.figure()
	plt.pcolor(X, Y, df, cmap='viridis')
	# plt.pcolor(df)
	plt.colorbar()
	# plt.clim(0, 30)
	plt.title(str(i))
	# plt.show()
	plt.savefig('data/plot'+str(i)+'.png')
	plt.close()

if os.path.isfile('out.mp4'):

	os.remove('out.mp4')

os.system('ffmpeg -f image2 -r 10 -i ' + datadir + 'plot%d.png -b 600k out.mp4')

os.system('totem out.mp4')