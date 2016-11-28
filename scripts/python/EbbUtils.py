import sys
import os
import struct
import numpy as np
from functools import partial

def importEbbSolution(slnFile):
	print("Importing " + slnFile)
	slnF = open(slnFile, 'rb')

	slnMagic = struct.unpack(">I", slnF.read(4))[0]
	slnVer = struct.unpack(">I", slnF.read(4))[0]
	nCells = struct.unpack(">I", slnF.read(4))[0]

	U = np.zeros((nCells, 4))

	t = struct.unpack(">d", slnF.read(8))[0]
	dt = struct.unpack(">d", slnF.read(8))[0]

	for idx in range(nCells):
		u1 = struct.unpack(">d", slnF.read(8))[0]
		u2 = struct.unpack(">d", slnF.read(8))[0]
		u3 = struct.unpack(">d", slnF.read(8))[0]
		u4 = struct.unpack(">d", slnF.read(8))[0]
		U[idx,:] = [u1, u2, u3, u4]

	return U

def importEbbMatlabMesh(meshFile):
	
	print("Importing " + meshFile)
	meshF = open(meshFile, 'rb')

	nNodes = struct.unpack(">Q", meshF.read(8))[0]
	V = np.zeros((nNodes, 2))
	print(str(nNodes) + " nodes")
	for idx in range(nNodes):
		x = struct.unpack(">d", meshF.read(8))[0]
		y = struct.unpack(">d", meshF.read(8))[0]
		V[idx,:] = [x, y]
	
	nEls = struct.unpack(">Q", meshF.read(8))[0]
	E = np.zeros((nEls, 3))
	print(str(nEls) + " elements")
	for idx in range(nEls):
		n1 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n2 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n3 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		E[idx,:] = [n1, n2, n3]

	ieSize = struct.unpack(">Q", meshF.read(8))[0]
	IE = np.zeros((ieSize, 4))
	print(str(ieSize) + " interior elements")
	for idx in range(ieSize):
		n1 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n2 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n3 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n4 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		IE[idx,:] = [n1, n2, n3, n4]

	beSize = struct.unpack(">Q", meshF.read(8))[0]
	BE = np.zeros((beSize, 4))
	print(str(beSize) + " boundary elements")
	for idx in range(beSize):
		n1 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n2 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n3 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		n4 = int(struct.unpack(">d", meshF.read(8))[0]) - 1
		BE[idx,:] = [n1, n2, n3, n4]

	mesh = {}
	mesh['V'] = V
	mesh['E'] = E
	mesh['BE'] = BE

	return mesh
	
def importEbbForces(forceFile):

	print("Importing " + forceFile)
	forcesF = open(forceFile, 'rb')
	fileSize = os.path.getsize(forceFile)
	
	numEntries = fileSize/(3*8)

	t = np.zeros((numEntries, 1))
	forces = np.zeros((numEntries, 2))
	i = 0
	with open(forceFile, 'rb') as forcesF:
		for data in iter(partial(forcesF.read, 3*8), ''):
			if len(data) == 3*8:
				t[i] = struct.unpack(">d", data[0:8])[0]
				forces[i,:] = [struct.unpack(">d", data[8:(8+8)])[0], struct.unpack(">d", data[(8+8):(8+8+8)])[0]]
				i = i+1
			else:
				return (t, forces)
