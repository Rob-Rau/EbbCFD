import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------
def getField(U, field):
	r = U[:,0]
	ru = U[:,1]
	rv = U[:,2]
	rE = U[:,3]

	g = 1.4
	s = field.lower()

	if (s == 'mach'):
		V = np.sqrt(ru**2 + rv**2)/r
		p = (g-1.)*(rE-0.5*r*V**2)
		c = np.sqrt(g*p/r)
		return V/c
	elif (s == 'pressure'):
		V = np.sqrt(ru**2 + rv**2)/r
		p = (g-1.)*(rE-0.5*r*V**2)
		return p
	elif (s == 'density'):
		return r
	elif (s == 'xmomentum'):
		return ru
	elif (s == 'ymomentum'):
		return rv
	elif (s == 'energy'):
		return rE/r
	elif (s == 'xvelocity'):
		return ru/r
	elif (s == 'yvelocity'):
		return rv/r


#-----------------------------------------------------------
def plotstate(Mesh, U, field, fname):
	V = Mesh['V']
	E = Mesh['E']
	BE = Mesh['BE']

	f = plt.figure(figsize=(12,6))

	F = getField(U, field)
	
	plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=F, shading='flat')

	for i in range(len(BE)):
		x = [V[BE[i,0],0], V[BE[i,1],0]]
		y = [V[BE[i,0],1], V[BE[i,1],1]]
		plt.plot(x, y, '-', linewidth=2, color='black')
	
	dosave = (len(fname) != 0)

	plt.axis('equal')

	plt.axis([-100, 100,-100, 100])
	#plt.axis([-2, 10,-4, 4])
	plt.colorbar()
	#plt.clim(0, 0.7)
	#plt.clim(9, 12)
	plt.title(field, fontsize=16)

	f.tight_layout()
	plt.show(block=(not dosave))
	if (dosave):
		plt.savefig(fname)
	
	plt.close(f)
