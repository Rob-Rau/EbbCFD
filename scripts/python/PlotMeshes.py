#!/usr/bin/env python3.4

import sys
import os

import matplotlib.pyplot as plt
import numpy as np

import EbbUtils as eu

#-----------------------------------------------------------
def plotmesh(Mesh, color):
	V = Mesh['V']
	E = Mesh['E']
	BE = Mesh['BE']

	F = np.zeros((np.size(E)/3, 1))
	F = F[:,0]
	F[:] = 1
	plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=F, edgecolors=color, cmap=plt.cm.gray, vmin=0, vmax=1, alpha=1, linewidth=1.4)

	for i in range(len(BE)):
		x = [V[BE[i,0],0], V[BE[i,1],0]]
		y = [V[BE[i,0],1], V[BE[i,1],1]]
		plt.plot(x, y, '-', linewidth=2, color=color)

	plt.axis('equal')
	plt.axis([-100, 100,-100, 100])


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print('Not enough input arguments')
		sys.exit()
	
	colors = ['k', 'b', 'g', 'c', 'y', 'm', 'r']

	f = plt.figure(figsize=(12,6))
	for idx in range(len(sys.argv) - 1):
		mesh = eu.importEbbMatlabMesh(sys.argv[idx+1])
		plotmesh(mesh, colors[idx%len(colors)])
		plt.hold(True)

	plt.show(block=True)