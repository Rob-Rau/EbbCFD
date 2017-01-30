#!/usr/bin/env python3

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
	
	if(F.shape[0] == V.shape[0]):
		plt.tripcolor(V[:,0], V[:,1], F, triangles=E, shading='gouraud', edgecolors=color, vmin=clim1, vmax=clim2, linewidth=1)
	else:
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

	meshIndicies = [i for i, s in enumerate(sys.argv) if(eu.meshTypeSupported(s))]
	slnIndicies = [i for i, s in enumerate(sys.argv) if 'sln' in s]

	meshFiles = []
	for idx in meshIndicies:
		meshFiles.append(sys.argv[idx])
	
	slnFiles = []
	for idx in slnIndicies:
		slnFiles.append(sys.argv[idx])

	sort_nicely(meshFiles)
	sort_nicely(slnFiles)

	globalMin = float("inf")
	globalMax = float("-inf")

	combineSln = False
	mesh = {}
	
	if((len(slnFiles) != len(meshFiles)) and (len(meshFiles) == 1)):
		combineSln = True
		mesh = eu.importMesh(meshFiles[0])
	elif(len(slnFiles) != len(meshFiles)):
		raise Exception("different number of solution files and mesh files! aborting")
		

	slns = []
	combinedSln = {}
	U = []
	edgeVals = []
	edges = []

	M = 0
	uAlloced = False
	for slnFile in slnFiles:
		sln = eu.importEbbSolution(slnFile)
		if((sln['version'] == 1) and combineSln):
			raise Exception("Cannot recombine solution files with file version 1")

		if(not combineSln):
			F = pu.getField(sln['U'], state)
			minVal = np.min(F)
			globalMin = min(globalMin, minVal)

			maxVal = np.max(F)
			globalMax = max(globalMax, maxVal)

			slns.append(sln)
		else:
			raise Exception("Recombing solutions not yet supported")

	print("globalMin = "+str(globalMin))
	print("globalMax = "+str(globalMax))
	
	globalMin = globalMin - 0.001
	globalMax = globalMax + 0.001

	colors = ['k', 'b', 'g', 'c', 'y', 'm', 'r']

	f = plt.figure(figsize=(12,6))

	for idx in range(len(meshFiles)):
		if(not combineSln):
			mesh = eu.importMesh(meshFiles[idx])
		print(slns[idx]['U'].shape)
		print(mesh['V'].shape)
		plotstate(mesh, slns[idx]['U'], state, "", globalMin, globalMax, colors[idx%len(colors)])
		plt.hold(True)

	plt.colorbar()
	plt.show(block=True)
	plt.close(f)
