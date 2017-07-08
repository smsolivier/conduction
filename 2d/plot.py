#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
import os 

datadir = 'data/'
prefix = 'out'

allfiles = os.listdir(datadir)

files = [] 

total = 0 

fmms = lambda x, y: np.sin(np.pi*x/.1) * np.sin(np.pi*y/.1)

for i in range(len(allfiles)):

	if (allfiles[i].startswith(prefix)):

		num = int(allfiles[i][len(prefix):])

		if (num > total): 

			total = num

N = total + 1
for i in range(N):

	print(i/N, end='\r')

	X = np.loadtxt('data/X'+str(i))
	Y = np.loadtxt('data/Y'+str(i))
	df = np.loadtxt('data/out'+str(i))

	plt.figure()
	# plt.pcolor(X, Y, df, cmap='viridis')
	plt.pcolor(X, Y, np.fabs(df - fmms(X, Y)), norm=LogNorm(), cmap='viridis')
	plt.colorbar()
	plt.title(str(i))
	plt.show()
	# plt.savefig('data/plot'+str(i)+'.png')
	# plt.close()

# if os.path.isfile('out.mp4'):

	# os.remove('out.mp4')

# os.system('ffmpeg -f image2 -r 10 -i ' + datadir + 'plot%d.png -b 600k out.mp4')

# os.system('totem out.mp4')