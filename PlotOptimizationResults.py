#!/usr/bin/env python3.4

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def PlotRealNetwork():
	weights = np.genfromtxt('SQPpoints_realNetwork.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0111292, 0.00998556, 0.00997867, 0.00972602, 0.00972354444459563524, 0.0141187, 0.0112997, 0.0112002, 0.010635, 0.0104412, 0.01034025185163719301, 0.0104114, 0.0105021, 0.010419, 0.011412, 0.0112074, 0.0108091]
	#12.6%

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	#plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')

	plt.suptitle('6 Processor Ethernet')
	plt.plot(f)

def PlotRealNetwork2():
	weights = np.genfromtxt('SQPpoints_realNetwork2.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0113527, 0.0100857, 0.00999296, 0.00991069, 0.00976466, 0.00974761, 0.00975821, 0.0138991, 0.0104268, 0.00999844, 0.00984646, 0.0107659, 0.0138946, 0.013898, 0.0108843, 0.0108483, 0.0115359, 0.0103734, 0.00991257, 0.0100628, 0.00978714814840781454, 0.0098085, 0.00980577, 0.00981618, 0.0110233, 0.00987941, 0.01355882222220922487, 0.010585, 0.00993282, 0.0098202]
	#13.5%

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	#plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')

	plt.suptitle('6 Processor Busy Ethernet')
	plt.plot(f)

def PlotNetworChange():
	weights = np.genfromtxt('SQPpoints_networkChange.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0180247, 0.0157878, 0.0145288, 0.0140524, 0.01389060726872196928, 0.0139973, 0.0139889, 0.0221232, 0.0178224, 0.0158195, 0.014936, 0.0141553, 0.0136331, 0.01335387053313078637]
	#22.4%

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	#plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')

	plt.suptitle('Simulated network change')
	plt.plot(f)

def PlotNetworChange2():
	weights = np.genfromtxt('SQPpoints_networkChange_p6.csv',
							skip_header=0,
							skip_footer=0,
							delimiter=',',
							dtype='float64')

	f = [0.0155781, 0.0147208, 0.0122766, 0.0119979, 0.0116224, 0.0114836, 0.0113497, 0.011347, 0.0184374, 0.0184385, 0.0174151, 0.0152859, 0.0141953, 0.013132, 0.0127488, 0.0121687, 0.0120014, 0.011813, 0.0116682, 0.0116633, 0.0113612, 0.01135814808033130778, 0.0113867, 0.0113881, 0.0113843]
	# 27.16%

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	#plt.xlabel('Optimizer iterations')
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
	# 25.422%
	
	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	#plt.xlabel('Optimizer iterations')
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
	# 10.789%

	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	#plt.xlabel('Optimizer iterations')
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

	# 29.28% speedup
	f = [0.0170872, 0.0149235, 0.0124973, 0.0122598, 0.01208770716631854016, 0.012091, 0.0120846]
	
	plt.subplot(2,1,1)
	#print(weights)
	plt.plot(weights)
	plt.title('Metis mesh weights')
	#plt.xlabel('Optimizer iterations')
	plt.ylabel('Mesh weights')

	plt.subplot(2,1,2)
	plt.title('Average iteration time')
	plt.xlabel('Optimizer iterations')
	plt.ylabel('Average solver iteration time [s]')
	plt.plot(f)

def PlotEbbSpeedup(elements, p, serialTime, parallelTime):
	speedup = serialTime/parallelTime
	handle = plt.plot(p, speedup, label=str(elements)+' element mesh')
	plt.axis('equal')
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

	matplotlib.rcParams.update({'font.size': 8})

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
	plt.savefig("Writeup/EbbSpeedup.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	h1 = PlotEbbEfficiency(6413, p, serialTime1, parallelTime1)
	h2 = PlotEbbEfficiency(11594, p, serialTime2, parallelTime2)
	plt.legend()
	plt.savefig("Writeup/EbbEfficiency.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	Plot6PevenSlowdown()
	plt.savefig("Writeup/6Pevenslow.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	Plot6PunevenSlowdown()
	plt.savefig("Writeup/6Punevenslow.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	Plot6PunevenSlowdownBigMesh()
	plt.savefig("Writeup/6PunevenslowBig.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	PlotNetworChange()
	plt.savefig("Writeup/5PnetworkChange.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	PlotRealNetwork()
	plt.savefig("Writeup/RealNetwork.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	PlotRealNetwork2()
	plt.savefig("Writeup/RealNetwork2.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.figure()
	PlotNetworChange2()
	plt.savefig("Writeup/6PNetworkChange.pdf", bbox_inches='tight', papertype='letter', format='pdf')

	plt.show(block=True)


