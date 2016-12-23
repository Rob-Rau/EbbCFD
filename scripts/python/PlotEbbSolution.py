#!/usr/bin/env python3.4

import re
import struct
import sys

import numpy as np
import matplotlib.pyplot as plt

import EbbUtils as eu
import PlotUtils as pu

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

#-----------------------------------------------------------
def plotstate(Mesh, U, field, fname, clim1, clim2, color):
	V = Mesh['V']
	E = Mesh['E']
	BE = Mesh['BE']

	#f = plt.figure(figsize=(12,6))
	F = pu.getField(U, field)
	
	#plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=F, edgecolors=color, shading='flat', vmin=clim1, vmax=clim2, linewidth=1)
	plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=F, shading='flat', vmin=clim1, vmax=clim2, linewidth=1)

	for i in range(len(BE)):
		x = [V[BE[i,0],0], V[BE[i,1],0]]
		y = [V[BE[i,0],1], V[BE[i,1],1]]
		plt.plot(x, y, '-', linewidth=2, color='black')
	
	dosave = (len(fname) != 0)

	plt.axis('equal')

	#plt.axis([-100, 100,-100, 100])
	#plt.axis([-2, 10,-4, 4])
	#plt.colorbar()
	#plt.clim(0, 0.7)
	#plt.clim(9, 12)
	plt.title(field, fontsize=16)

	#f.tight_layout()
	#plt.show(block=(not dosave))
	#if (dosave):
	#	plt.savefig(fname)
	
	#plt.close(f)

#-----------------------------------------------------------
if __name__ == "__main__":
	if len(sys.argv) < 4:
		print('Not enough input arguments')
		sys.exit()

	state = sys.argv[1]

	meshIndicies = [i for i, s in enumerate(sys.argv) if 'mmsh' in s]
	slnIndicies = [i for i, s in enumerate(sys.argv) if 'sln' in s]

	meshFiles = []
	for idx in meshIndicies:
		meshFiles.append(sys.argv[idx])
	
	slnFiles = []
	for idx in slnIndicies:
		slnFiles.append(sys.argv[idx])

	sort_nicely(meshFiles)
	sort_nicely(slnFiles)

#	meshFile = sys.argv[1]
#	slnFile = sys.argv[2]
#	state = sys.argv[3]
#	plotSave = ""
#
#	if len(sys.argv) == 5:
#		plotSave = sys.argv[4]

#	mesh = eu.importEbbMatlabMesh(meshFile)
#	U = eu.importEbbSolution(slnFile)
#
#	pu.plotstate(mesh, U, state, plotSave)

	globalMin = float("inf")
	globalMax = float("-inf")

	slns = []
	for slnFile in slnFiles:
		sln = eu.importEbbSolution(slnFile)
		F = pu.getField(sln, state)
		minVal = np.min(F)
		print("minVal = "+str(minVal))
		globalMin = min(globalMin, minVal)

		maxVal = np.max(F)
		print("maxVal = "+str(maxVal))
		globalMax = max(globalMax, maxVal)

		slns.append(sln)

	globalMin = globalMin - 0.001
	globalMax = globalMax + 0.001

	print("globalMin = "+str(globalMin))
	print("globalMax = "+str(globalMax))
	colors = ['k', 'b', 'g', 'c', 'y', 'm', 'r']

	f = plt.figure(figsize=(12,6))
	
	# globalMin = 0.45
	# globalMax = 1.8

	for idx in range(len(meshFiles)):
		mesh = eu.importEbbMatlabMesh(meshFiles[idx])
		#sln = eu.importEbbSolution(slnFiles[idx])
		plotstate(mesh, slns[idx], state, "", globalMin, globalMax, colors[idx%len(colors)])
		plt.hold(True)

	#plt.axis([-0.3, 1.3,-0.5, 0.5])

	plt.colorbar()
	plt.show(block=True)
	plt.close(f)