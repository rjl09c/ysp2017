import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder
from pylab import*
from numpy import ma

#derives velx with respect to x
def derivx(vel,xcoords):

	x1 = 0
	x2 = 0
	x3 = 0
	distance = xcoords[1][0] - xcoords[0][0]
	velxdx = np.zeros((320,320))
	for i in range(len(vel)):
		for x in range(len(vel)):
			if 0 < i < len(vel) - 1:
				velxdx[i,x] = ((-1/2) * vel[i - 1][x]) + ((1/2) * vel[i + 1][x])
			elif i == 0:
				velxdx[i,x] = (((-3/2) * vel[i][x]) + (2 * vel[i + 1][x]) + ((-1/2) * vel[i+2][x]))
			elif i == len(vel) - 1:
				velxdx[i,x] = ((-3/2) * vel[i][x]) + (2 * vel[i - 1][x]) + ((-1/2) * vel[i-2][x])

	return velxdx/distance


#derives vely with respect to y
def derivy(vel,xcoords):

	y1 = 0
	y2 = 0
	y3 = 0
	distance = xcoords[1][0] - xcoords[0][0]
	velydy = np.zeros((320,320))
	for i in range(len(vel)):
		for x in range(len(vel)):
			if 0 < x < len(vel) - 1:
				velydy[i,x] = (((-1/2) * vel[i][x - 1]) + ((1/2) * vel[i][x + 1]))
			elif x == 0:
				velydy[i,x] = (((-3/2)*vel[i][x]) + (2*vel[i][x+1])+((-1/2) * vel[i][x + 2]))
			elif x == len(vel) - 1:
				velydy[i,x] = (((-3/2)*vel[i][x]) + (2*vel[i][x - 1]) + ((-1/2) * vel[i][x - 2]))

	return velydy/distance

#find magnitude
def magnitude(velx,vely, xcoords):
	mag = np.zeros((320,320))
	yderiv = derivy(vely, xcoords)
	xderiv = derivx(velx, xcoords)
	for i in range(len(xderiv)):
		for x in range(len(xderiv)):
			mag[i][x] = log((((yderiv[i,x]**2) + (xderiv[i,x]**2))**.5))

	return mag

#quiver style plots
def quiver(xcoords,ycoords,velx,vely):
	X = xcoords
	Y = ycoords
	U = derivx(velx, xcoords)
	V = derivy(vely, xcoords)
	C = magnitude(velx, vely, xcoords)
	#R=magnitude(velx, vely, xcoords)
	'''
	plt.figure()
	plt.title('Arrows scale with plot width, not view')
	Q = plt.quiver(X, Y, U, V, units='width')
	qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
	                   coordinates='figure')
	'''
	plt.figure()
	#M = np.hypot(U, V)
	Q = plt.quiver(X, Y, U, V, units = 'xy', angles = 'xy', pivot = 'tip',
                   scale = 1/0.1, minlength = 0)
	qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
			           coordinates='figure')
	plt.scatter(X[::3, ::3], Y[::3, ::3], s=5)
	'''
	plt.figure()
	plt.title("Magnitude of the Velocity Gradient")
	M = np.hypot(U, V)
	Q = plt.quiver(xcoords, ycoords, U, V, M, units='x', pivot='tip', width=0.022,
	               scale=1 / 0.15)
	#qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
	#                   coordinates='figure')

	plt.scatter(xcoords, ycoords, c=magnitude(velx, vely, xcoords), s=5)
	'''
	#plt.gcf().set_size_inches(60,50)
	plt.show()



ds = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0000")
ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge, 
									dims=ds.domain_dimensions)

xcoords = np.array(ad["x"])
ycoords = np.array(ad["y"])
velx = np.array(ad["velx"])
vely = np.array(ad["vely"])


ds1 = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0001")
dd = ds1.covering_grid(level=0, left_edge=ds1.index.grids[0].LeftEdge, dims=ds1.domain_dimensions)

xcoords1 = np.array(dd["x"])
ycoords1 = np.array(dd["y"])
velx1 = np.array(dd["velx"])
vely1 = np.array(dd["vely"])

'''
plt.figure()
plt.scatter(xcoords, ycoords,c=magnitude(velx, vely, xcoords), marker= 'o',edgecolor='none')
cb = plt.colorbar()
cb.set_label('Log of the Magnitude')

plt.figure
plt.scatter(xcoords1, ycoords1,c=finalvorticity(velx1,vely1, xcoords1), marker= 'o',edgecolor='none')
cb = plt.colorbar()
cb.set_label('Log of the Magnitude')

plt.show()
'''
quiver(xcoords1,ycoords1,velx1,vely1)
