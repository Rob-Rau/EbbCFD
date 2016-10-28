#!/usr/bin/env python3.4

import struct
import sys

import EbbUtils as eu
import PlotUtils as pu

#-----------------------------------------------------------
if __name__ == "__main__":
	if len(sys.argv) < 4:
		print('Not enough input arguments')
		sys.exit()

	meshFile = sys.argv[1]
	slnFile = sys.argv[2]
	state = sys.argv[3]
	plotSave = ""

	if len(sys.argv) == 5:
		plotSave = sys.argv[4]

	mesh = eu.importEbbMatlabMesh(meshFile)
	U = eu.importEbbSolution(slnFile)

	pu.plotstate(mesh, U, state, plotSave)
