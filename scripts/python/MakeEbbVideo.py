#!/usr/bin/env python3.4

import re
import sys
import os
import subprocess
#from subprocess import run
from os import listdir
from os.path import isfile, join

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import EbbUtils as eu
import PlotUtils as pu

def plotstate(Mesh, U, field, fname):
	V = Mesh['V']
	E = Mesh['E']
	BE = Mesh['BE']

	f = plt.figure(figsize=(12,6))

	F = pu.getField(U, field)
	plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=F, shading='flat', vmin=0.0, vmax=0.8)

	for i in range(len(BE)):
		x = [V[BE[i,0],0], V[BE[i,1],0]]
		y = [V[BE[i,0],1], V[BE[i,1],1]]
		plt.plot(x, y, '-', linewidth=2, color='black')
	
	#dosave = (len(fname) != 0)

	plt.axis('equal')

	#plt.axis([-100, 100,-100, 100])
	plt.axis([-2, 10,-4, 4])
	plt.colorbar()
	#plt.clim(0, 0.8)
	plt.title(field, fontsize=16)

	f.tight_layout()
	plt.show()#block=(not dosave))
	#if (dosave):
	plt.savefig(fname)
	
	plt.close(f)


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

if __name__ == "__main__":
	
	framerate = 30
	if len(sys.argv) < 3:
		print('Not enough input arguments')
		sys.exit()

	meshFile = sys.argv[1]
	state = sys.argv[2]

	if len(sys.argv) < 4:
		print('Assuming 30 fps')
	elif len(sys.argv) == 4:
		framerate = int(sys.argv[3])

	mesh = eu.importEbbMatlabMesh(meshFile)

	slnFiles = [f for f in listdir("./") if (isfile(join("./", f)) and ("esln" in f) and ("final" not in f))]
	sort_nicely(slnFiles)
	
	os.mkdir("output")
	for idx in range(len(slnFiles)):
		U = eu.importEbbSolution(slnFiles[idx])
		plotstate(mesh, U, state, "output/"+str(idx).zfill(8)+".png")
	
	try:
		subprocess.call(["ffmpeg", "-framerate", str(framerate), "-i", "output/%08d.png", "-c:v", "libx264", "output/"+state+".mp4"], stdout=subprocess.PIPE)
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			print("ffmpeg not installed, trying avconv")
			subprocess.call(["avconv", "-framerate", str(framerate), "-i", "output/%08d.png", "-c:v", "libx264", "output/"+state+".mp4"], stdout=subprocess.PIPE)
		else:
			print("Failed to convert video")

