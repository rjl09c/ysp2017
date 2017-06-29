import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder
from pylab import*
from numpy import ma
from numpy import linalg as LA


#derives vel with respect to x
def derivx(vel,xcoords):
	distance = xcoords[1][0] - xcoords[0][0]
	velxdx = np.zeros((320,320))
	
	for i in range(len(vel)):
		
		for x in range(len(vel)):
			
			if 0 < i < len(vel) - 1:
				velxdx[i,x] = ((-1/2) * vel[i-1][x]) + ((1/2) * vel[i+1][x])
			
			elif i == 0:
				velxdx[i,x] = (((-3/2) * vel[i][x]) + (2 * vel[i+1][x]) + ((-1/2) * vel[i+2][x]))
			
			elif i == len(vel) - 1:
				velxdx[i,x] = ((-3/2) * vel[i][x]) + (2 * vel[i-1][x]) + ((-1/2) * vel[i-2][x])

	return velxdx/distance


#derives vel with respect to y
def derivy(vel,xcoords):
	distance = xcoords[1][0] - xcoords[0][0]
	velydy = np.zeros((320,320))

	for i in range(len(vel)):

		for x in range(len(vel)):

			if 0 < x < len(vel) - 1:
				velydy[i,x] = (((-1/2) * vel[i][x-1]) + ((1/2) * vel[i][x+1]))
			
			elif x == 0:
				velydy[i,x] = (((-3/2)*vel[i][x]) + (2*vel[i][x+1]) + ((-1/2) * vel[i][x + 2]))
			
			elif x == len(vel) - 1:
				velydy[i,x] = (((-3/2)*vel[i][x]) + (2*vel[i][x-1]) + ((-1/2) * vel[i][x-2]))

	return velydy/distance


#second derivative of vel with respect to x
def deriv2x(vel,xcoords):
	distance = xcoords[1][0] - xcoords[0][0]
	velxdx = np.zeros((320,320))

	for i in range(len(vel)):

		for x in range(len(vel)):

			if 0 < i < len(vel) - 1:
				velxdx[i,x] = (vel[i-1][x]) + (-2 * vel[i][x]) + (vel[i+1][x])

			elif i == 0:
				velxdx[i,x] = ((2 * vel[i][x]) + (-5 * vel[i+1][x]) + (4* vel[i+2][x]) + (-1 * vel[i+3][x]))

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
				mixed[i][x] = ((vel[i+1][x+1]) - (vel[i+1][x-1]) - (vel[i-1][x+1]) + (vel[i-1][x-1]))/(4*distx*disty)

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
	allhessian = [[[] for j in range(320)] for i in range(320)]
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

				allhessian[j][k] = hessian
	allhessian = np.array(allhessian)
	return allhessian


#find determinant
def determinant(allhessian):
	deters = np.zeros((320,320))

	for j in range(len(allhessian)):

		for k in range(len(allhessian)):

			x = allhessian[j,k]

			deters[j,k] = (x[0,0]*x[1,1]) - (x[1,0]*x[0,1])

	return deters


#find magnitude
def magnitude(velx,vely, xcoords):
	mag = np.zeros((320,320))
	yderiv = derivy(vely, xcoords)
	xderiv = derivx(velx, xcoords)

	for i in range(len(xderiv)):

		for x in range(len(xderiv)):
			mag[i][x] = (((yderiv[i,x]**2) + (xderiv[i,x]**2))**.5)

	return mag


#finds extrema and saddlepoints
def extrema(allhessian, velx, vely, xcoords):
	deters = determinant(allhessian)
	extrem = np.zeros((320,320))
	mag = magnitude(velx, vely, xcoords)

	for j in range(len(extrem)):

		for k in range(len(extrem)):

			if mag[j][k] == 0:
				if deters[j,k] < 0:
					extrem[j, k] = -1

				elif deters[j,k] == 0:
					extrem[j,k] = 0

				else:
					x = allhessian[j,k]
					if deter[j,k] > 0 and x[0,0] > 0:
						extem[j, k] = -2
					elif deter[j,k] > 0 and x[0,0] < 0:
						extrem[j, k] = 2
	return extrem				


#creates jacobia matrix for each point						
def jacobian(xcoords,velx, vely):
	xx = derivx(velx, xcoords)
	xy = derivy(velx, xcoords)
	yx = derivx(vely, xcoords)
	yy = derivy(vely, xcoords)
	jacob = np.zeros ((2,2))
	alljacob = [[[] for j in range(320)] for i in range(320)]

	for j in range(len(xx)):

		for k in range(len(xx)):

			for i in range(len(jacob)):

				for c in range(len(jacob)):

					if c == 0 and i == 0:
						jacob [i][c] = xx[j][k] 

					elif c == 1 and i == 0:
						jacob[i][c] = xy[j][k]

					elif c ==1 and i == 1:
						jacob[i][c] = yy[j][k]

				alljacob[j][k] = jacob
	alljacob = np.array(alljacob)			
	return alljacob


#obtains eigenvalues for all points' jacobian matrices and then checks the extrema
def evals(alljacob):
	eigen = [[[] for j in range(320)] for i in range(320)]
	extrema = np.zeros((320,320))

	for j in range(len(alljacob)):

		for k in range(len(alljacob)):

			x = alljacob[j,k]
			eigen[j][k] = LA.eigvalsh(x)
			y = eigen [j][k]

			if y[0]>0 and y[1] > 0:
				extrema[j,k] = 2

			elif y[0]<0 and y[1] <0:
				extrema[j,k] = -2

			elif y[0]*y[1] < 0:
				extrema[j,k] = 3

	return extrema

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
	extrema(hess(xcoords, ycoords, velx), velx, vely, xcoords)

	#creates Hessian marix for y velocity for file 1
	(extrema(hess(xcoords, ycoords, vely), velx, vely, xcoords))

	#prints extrema for file1
	print(evals(jacobian(xcoords, velx, vely)))
	'''plt.figure()
	plt.scatter(xcoords, ycoords,c=evals(jacobian(xcoords, velx, vely)), marker= 'o',edgecolor='none')
	cb = plt.colorbar()
	cb.set_label('Extrema')
	plt.show()'''


main()
