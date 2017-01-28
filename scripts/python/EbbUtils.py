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
	t = struct.unpack(">d", slnF.read(8))[0]
	dt = struct.unpack(">d", slnF.read(8))[0]

	order = 1
	sln = {}
	sln['version'] = slnVer

	if(slnVer == 1):
		sln['M'] = 4
		print("version 1")
		U = np.zeros((nCells, 4))

		for idx in range(nCells):
			u = []
			for ui in range(4):
				u.append(struct.unpack(">d", slnF.read(8))[0])
			
			U[idx,:] = u

		sln['U'] = U

	elif(slnVer == 2):
		print("version 2")
		M = struct.unpack(">I", slnF.read(4))[0]
		sln['M'] = M
		order = struct.unpack(">I", slnF.read(4))[0]

		U = np.zeros((nCells, M))
		gens = [] # global element numbers
		edges = []
		edgeVals = []

		for idx in range(nCells):
			# read global cell number and number of faces
			gens.append(struct.unpack(">I", slnF.read(4))[0])
			edges.append(struct.unpack(">B", slnF.read(1))[0])
			# read cell average value
			u = []
			for ui in range(M):
				u.append(struct.unpack(">d", slnF.read(8))[0])
			
			U[idx,:] = u

			edgeVal = []
			#read edge values
			for en in range(edges[-1]):
				u = []
				for ui in range(M):
					u.append(struct.unpack(">d", slnF.read(8))[0])

				edgeVal.append(u)

			edgeVals.append(edgeVal)

		sln['U'] = U
		sln['gens'] = gens
		sln['edges'] = edges
		sln['edgeVals'] = edgeVals

	else:
		raise Exception("Solution file version "+str(slnVer)+" is not supported")

	sln['order'] = order
	return sln

def meshTypeSupported(meshFile):
	if("mmsh" in meshFile):
		return True
	else:
		return False

def importMesh(meshFile):
	if("mmsh" in meshFile):
		return importEbbMatlabMesh(meshFile)
	else:
		raise Exception("Cannot import "+meshFile+". Unsupported file type")

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
	mesh['IE'] = IE

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
