import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder
from pylab import*

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

ds = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0000")
ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge, 
									dims=ds.domain_dimensions)

xcoords = ad["x"]
ycoords = ad["y"]
velx = ad["velx"]
vely = ad["vely"]


ds1 = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0001")
dd = ds1.covering_grid(level=0, left_edge=ds1.index.grids[0].LeftEdge, dims=ds1.domain_dimensions)

xcoords1 = dd["x"]
ycoords1 = dd["y"]
velx1 = dd["velx"]
vely1 = dd["vely"]


plt.scatter(xcoords, ycoords,c=magnitude(velx, vely, xcoords), marker= 'o',edgecolor='none')
cb = plt.colorbar()
cb.set_label('Magnitude')
plt.show()
