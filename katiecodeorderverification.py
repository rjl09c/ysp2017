import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder
from pylab import*
from numpy import ma
from numpy import linalg as LA


#deriveswith respect to x
def derivx(vel,xcoords):
	distance = xcoords[1][0] + xcoords[0][1] - 2*xcoords[0][0]
	velxdx = np.zeros((200,200))
	
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
	distance = xcoords[1][0] + xcoords[0][1] - 2*xcoords[0][0]
	velydy = np.zeros((200,200))

	for i in range(len(vel)):

		for x in range(len(vel)):

			if 0 < x < len(vel) - 1:
				velydy[i,x] = (((-1/2) * vel[i][x-1]) + ((1/2) * vel[i][x+1]))
			
			elif x == 0:
				velydy[i,x] = (((-3/2)*vel[i][x]) + (2*vel[i][x+1]) + ((-1/2) * vel[i][x + 2]))
			
			elif x == len(vel) - 1:
				velydy[i,x] = (((-3/2)*vel[i][x]) + (2*vel[i][x-1]) + ((-1/2) * vel[i][x-2]))

	return velydy/distance

#calculating l2 norm for csv files
def  normFile1(ycoords, velx, velx1):
	e1 = 0

	norm = np.zeros((100,100))

	for i in range(len(ycoords)):

		for j in range(len(ycoords)):

			norm [i][j] = (abs(float(velx[i][j])-float(velx1[i][j])))
			e1 = e1 + (abs((velx[i][j])-float(velx1[i][j])))

	e1 = (e1/len(velx1))
	#return norm
	return e1


#calculating l2 norm for csv files
def  normFile2(ycoords, velx, velx1):
	e2 = 0

	norm = np.zeros((200,200))

	for i in range(len(ycoords)):

		for j in range(len(ycoords)):

			norm [i][j] = (abs(float(velx[i][j])-float(velx1[i][j])))
			e2 = e2 + (abs(float(velx[i][j])-float(velx1[i][j])))

	e2 = (e2/len(velx1))
	#return norm
	return e2


#second derivative of vel with respect to x
def deriv2x(vel,xcoords):
	distance = xcoords[1][0] - xcoords[0][0]
	velxdx = np.zeros((100,100))

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
	velydy = np.zeros((100,100))

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
	mixed = np.zeros((100,100))
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
	allhessian = [[[] for j in range(100)] for i in range(100)]
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
	deters = np.zeros((100,100))

	for j in range(len(allhessian)):

		for k in range(len(allhessian)):

			x = allhessian[j,k]

			deters[j,k] = (x[0,0]*x[1,1]) - (x[1,0]*x[0,1])

	return deters


#find magnitude
def magnitude(velx,vely, xcoords):
	mag = np.zeros((100,100))
	yderiv = derivy(vely, xcoords)
	xderiv = derivx(velx, xcoords)

	for i in range(len(xderiv)):

		for x in range(len(xderiv)):
			mag[i][x] = (((yderiv[i,x]**2) + (xderiv[i,x]**2))**.5)

	return mag


#finds extrema and saddlepoints
def extrema(allhessian, velx, vely, xcoords):
	deters = determinant(allhessian)
	extrem = np.zeros((100,100))
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
	alljacob = [[[] for j in range(100)] for i in range(100)]

	for j in range(len(alljacob)):

		for k in range(len(alljacob)):

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
	eigen = [[[] for j in range(100)] for i in range(100)]
	extrema = np.zeros((100,100))

	for j in range(len(alljacob)):

		for k in range(len(alljacob)):

			x = alljacob[j,k]
			eigen[j][k] = LA.eigvalsh(x)
			y = eigen [j][k]

			if y[0]>0 and y[1]>0:
				extrema[j,k] = 2

			elif y[0]<0 and y[1]<0:
				extrema[j,k] = -2

			elif y[0]*y[1]<0:
				extrema[j,k] = 3

	return extrema

#main function
def main():

	zvals1 = np.loadtxt("Grid1.csv", dtype ='float', delimiter = ',')
	zvals2 = np.loadtxt("Grid2.csv", dtype ='float', delimiter = ',')
	xyvals1 = np.loadtxt("xy1.csv", dtype ='float', delimiter = ',')
	xyvals2 = np.loadtxt("xy2.csv", dtype ='float', delimiter = ',')
	dx1vals = np.loadtxt("Dx1.csv", dtype ='float', delimiter = ',')
	dy1vals = np.loadtxt("Dy1.csv", dtype ='float', delimiter = ',')
	dx2vals = np.loadtxt("Dx2.csv", dtype ='float', delimiter = ',')
	dy2vals = np.loadtxt("Dy2.csv", dtype ='float', delimiter = ',')
	
	#norms grid1
	x = np.meshgrid(xyvals1, xyvals1)[0]
	y = np.meshgrid(xyvals1, xyvals1)[1]
	normnewfilex = normFile1(y, dx1vals, derivx(zvals1, x))
	normnewfiley = normFile1(y, dy1vals, derivy(zvals1, x))

	#norms grid2
	x1 = np.meshgrid(xyvals2, xyvals2)[0]
	y1 = np.meshgrid(xyvals2, xyvals2)[1]
	normnewfilex1 = normFile2(y1, dx2vals, derivx(zvals2, x1))
	normnewfiley1 = normFile2(y1, dy2vals, derivy(zvals2, x1))

	#graphs of norms

	#grid1 norms 
	'''
	plt.figure()

	plt.scatter(x, y, c = normnewfilex, marker= 'o',edgecolor='none')
	plt.colorbar()
	plt.show()

	plt.figure()

	plt.scatter(x, y, c = normnewfiley, marker= 'o',edgecolor='none')
	plt.colorbar()
	plt.show()
	'''

	#grid2 norms
	'''
	plt.figure()

	plt.scatter(x1, y1, c = normnewfilex1, marker= 'o',edgecolor='none')
	plt.colorbar()
	plt.show()

	plt.figure()

	plt.scatter(x1, y1, c = normnewfiley1, marker= 'o',edgecolor='none')
	plt.colorbar()
	plt.show()
	'''

	#dx error norms as a function of h
	dxnorm = log(normnewfilex/normnewfilex1)/log((xyvals1[1]-xyvals1[0])/(xyvals2[1]-xyvals2[0]))
	print(abs(dxnorm))

	#dy error norms as a function of h 
	dynorm = log(normnewfiley/normnewfiley1)/log((xyvals1[1]-xyvals1[0])/(xyvals2[1]-xyvals2[0]))
	print(abs(dynorm))
	
	#prints extrema for file1 using jacobian method
	print(evals(jacobian(x, zvals1, zvals1)))

	#prints extrema for file1 using hessian method and second derivatives (which are missing at the moment)
	#(extrema(hess(x, y, zvals1), zvals1, zvals1, x))


main()
