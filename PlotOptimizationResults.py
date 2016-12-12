#!/usr/bin/env python3.4

import sys
import numpy as np
import matplotlib.pyplot as plt

def PlotRealNetwork():
	weights = np.genfromtxt('SQPpoints_realNetwork.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0111292, 0.00998556, 0.00997867, 0.00972602, 0.00972354444459563524, 0.0141187, 0.0112997, 0.0112002, 0.010635, 0.0104412, 0.01034025185163719301, 0.0104114, 0.0105021, 0.010419, 0.011412, 0.0112074, 0.0108091]

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')
	plt.plot(f)

def PlotNetworChange():
	weights = np.genfromtxt('SQPpoints_networkChange.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0180247, 0.0157878, 0.0145288, 0.0140524, 0.01389060726872196928, 0.0139973, 0.0139889, 0.0221232, 0.0178224, 0.0158195, 0.014936, 0.0141553, 0.0136331, 0.01335387053313078637]

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')
	plt.plot(f)

def PlotNetworChange2():
	weights = np.genfromtxt('SQPpoints_networkChange_p6.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0155781, 0.0147208, 0.0122766, 0.0119979, 0.0116224, 0.0114836, 0.0113497, 0.011347, 0.0184374, 0.0184385, 0.0174151, 0.0152859, 0.0141953, 0.013132, 0.0127488, 0.0121687, 0.0120014, 0.011813, 0.0116682, 0.0116633, 0.0113612, 0.01135814808033130778, 0.0113867, 0.0113881, 0.0113843]

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')
	plt.plot(f)



def Plot6PevenSlowdown():
	weights = np.genfromtxt('SQPpoints_p6_evenSlow.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0174039, 0.01601, 0.0130469, 0.01292756310215702718, 0.0129665, 0.0129794]
	
	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')
	plt.plot(f)

def Plot6PunevenSlowdownBigMesh():
	weights = np.genfromtxt('SQPpoints_p6_unevenSlow_bigmesh.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0376879, 0.033933, 0.0339323, 0.0339995, 0.0339766, 0.033817, 0.033716, 0.0337312, 0.03370082908206516181, 0.0336107, 0.0336217]
	
	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')
	plt.plot(f)

def Plot6PunevenSlowdown():
	weights = np.genfromtxt('SQPpoints_p6_unevenSlow.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	# 29.26% speedup
	f = [0.0170872, 0.0149235, 0.0124973, 0.0122598, 0.01208770716631854016, 0.012091, 0.0120846]
	
	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')
	plt.plot(f)

def PlotEbbSpeedup(elements, p, serialTime, parallelTime):
	speedup = serialTime/parallelTime
	handle = plt.plot(p, speedup, label=str(elements)+' element mesh')
	plt.title('speedup')
	return handle

def PlotEbbEfficiency(elements, p, serialTime, parallelTime):
	efficiency = serialTime/(p*parallelTime)
	handle = plt.plot(p, efficiency, label=str(elements)+' element mesh')
	plt.title('efficiency')
	return handle

if __name__ == "__main__":
	
	p = np.zeros((5, 1))
	p[0] = 2
	p[1] = 4
	p[2] = 8
	p[3] = 16
	p[4] = 32

	parallelTime1 = np.zeros((5, 1))
	serialTime1 = 150.15
	parallelTime1[0] = 73.079
	parallelTime1[1] = 38.918
	parallelTime1[2] = 19.652
	parallelTime1[3] = 11.154
	parallelTime1[4] = 8.0195

	parallelTime2 = np.zeros((5, 1))
	serialTime2 = 272.58
	parallelTime2[0] = 138.07
	parallelTime2[1] = 65.637
	parallelTime2[2] = 35.804
	parallelTime2[3] = 19.787
	parallelTime2[4] = 12.877

	plt.figure()
	h1 = PlotEbbSpeedup(6413, p, serialTime1, parallelTime1)
	h2 = PlotEbbSpeedup(11594, p, serialTime2, parallelTime2)
	plt.legend()

	plt.figure()
	h1 = PlotEbbEfficiency(6413, p, serialTime1, parallelTime1)
	h2 = PlotEbbEfficiency(11594, p, serialTime2, parallelTime2)
	plt.legend()

	plt.figure()
	Plot6PevenSlowdown()

	plt.figure()
	Plot6PunevenSlowdown()

	plt.figure()
	Plot6PunevenSlowdownBigMesh()

	plt.figure()
	PlotNetworChange()

	plt.figure()
	PlotRealNetwork()

	plt.figure()
	PlotNetworChange2()
	
	plt.show(block=True)


