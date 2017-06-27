import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder
from pylab import*
from numpy import ma


#derives vel with respect to x
def derivx(vel,xcoords):
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


#derives vel with respect to y
def derivy(vel,xcoords):
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


#second derivative of vel with respect to x
def deriv2x(vel,xcoords):
	distance = xcoords[1][0] - xcoords[0][0]
	velxdx = np.zeros((320,320))

	for i in range(len(vel)):

		for x in range(len(vel)):

			if 0 < i < len(vel) - 1:
				velxdx[i,x] = (vel[i - 1][x]) + (-2 * vel[i][x]) + (vel[i + 1][x])

			elif i == 0:
				velxdx[i,x] = ((2 * vel[i][x]) + (-5 * vel[i + 1][x]) + (4* vel[i+2][x]) + (-1 * vel[i+3][x]))

			elif i == len(vel) - 1:
				velxdx[i,x] = ((-3/2) * vel[i][x]) + (2 * vel[i-1][x]) + ((-1/2) * vel[i-2][x])

	return velxdx/distance


#second derivative of vel with respect to y
def deriv2y(vel,xcoords):
	distance = xcoords[1][0] - xcoords[0][0]
	velydy = np.zeros((320,320))

	for i in range(len(vel)):

		for x in range(len(vel)):

			if 0 < x < len(vel) - 1:
				velydy[i,x] = ((vel[i][x-1]) + (-2 * vel[i][x]) + (vel[i][x+1]))

			elif x == 0:
				velydy[i,x] = (((2)*vel[i][x]) + (-5 * vel[i][x+1]) + ((4) * vel[i][x+2]) + (-1 * vel[i][x+3]))

			elif x == len(vel) - 1:
				velydy[i,x] = (((2) * vel[i][x]) + (-5 * vel[i][x - 1]) + ((4) * vel[i][x-2]) + (-1 * vel[i][x-3]))

	return velydy/distance


#second derivative of a mixed derivative
def mixed_deriv(xcoords, ycoords, vel):
	distx = xcoords[1][0] - xcoords[0][0]
	disty = ycoords[0][1] - ycoords[0][0]
	mixed = np.zeros((320,320))
	veldx = derivx(vel, xcoords)
	veldy = derivy(veldx, xcoords)  #takes deriv of vel with respect to x and derives that in the y direction

	for i in range(len(vel)):

		for x in range(len(vel)):

			if 0 < i < len(vel) - 1 and 0 < x < len(vel) - 1:
				mixed[i][x] = ((vel[i+1][x + 1]) - (vel[i+1][x-1]) - (vel[i - 1][x + 1]) + (vel[i -1][x-1]))/(4*distx*disty)

			#if on edges derives with respect to x first
			elif i == 0 or i == len(vel) - 1 or x == 0 or x == len(vel) - 1:
				mixed[i][x]=veldy[i][x]

	return mixed

#create hessian matrix for each point
def hess(xcoords, ycoords, vel):
	veldx = deriv2x(vel, xcoords)  #retrieves the second derivatives of the velocity in the x direction
	veldy = deriv2y(vel, xcoords)  #retrieves the second derivatives of the velocity in the y direction
	mixed = mixed_deriv(xcoords, ycoords, vel) #retrieves the second mixed derivatives of the velocity 
	hessian = np.zeros((2,2))

	for j in range(len(veldx)):

		for k in range(len(veldx)):

			for i in range(len(hessian)):

				for x in range(len(hessian)):

					if i == 0 and x == 1:
						hessian[i,x] = mixed[j,k]
						hessian[i+1][x-1] = mixed[j,k]

					elif x == 0 and i == 0:
						hessian[i,x] = veldx[j,k]

					elif x == 1 and i == 1:
						hessian[i,x] = veldy[j,k] 

			if j == 0 and k < 100:
				print(hessian)


#loads files and calls hess function
def main():
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

	#creates Hessian matrix for x velocity for file 1
	hess(xcoords, ycoords, velx)

	#creates Hessian marix for y velocity for file 1
	#hess(xcoords, ycoords, vely)


main()
