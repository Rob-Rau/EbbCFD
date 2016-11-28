#!/usr/bin/env python3.4

import sys
import EbbUtils as eu
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
	
	if len(sys.argv) < 3:
		print("Not enough input args, exiting")
		sys.exit()

	forceFile = sys.argv[1]
	forceComp = int(sys.argv[2])

	if (forceComp > 1) or (forceComp < 0):
		print("Force component should be 0 or 1")
		sys.exit()

	(t, forces) = eu.importEbbForces(forceFile)

	force = forces[:,forceComp]

	plt.plot(t, force)
	plt.show(block=True)
