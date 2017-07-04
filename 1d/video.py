#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import os 

datadir = 'data/'
prefix = 'out_'

allfiles = os.listdir(datadir)

files = [] 

total = 0 

for i in range(len(allfiles)):

	if (allfiles[i].startswith(prefix)):

		num = int(allfiles[i][len(prefix):])

		if (num > total): 

			total = num

for i in range(total+1):

	print(i/(total+1), end='\r')

	x, y = np.loadtxt(datadir+prefix+str(i), unpack=True)

	plt.figure()
	plt.plot(x, y)
	plt.ylim(0, 10)
	plt.title(str(i))
	plt.savefig(datadir+'plot'+str(i)+'.png')
	plt.close()

os.system('ffmpeg -f image2 -r 10 -i ' + datadir + 'plot%d.png -b 600k out.mp4') 