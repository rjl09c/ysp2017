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
def  norm(ycoords, velx, velx1):
	sumofy = 0

	norm = np.zeros((200,200))

	for i in range(len(ycoords)):

		for j in range(len(ycoords)):

			norm [i][j] = (abs(float(velx[i][j])-float(velx1[i][j])))
			sumofy = sumofy + ((float(velx[i][j])-float(velx1[i][j]))**2)**.5

	sumofy = (sumofy/len(velx1))
	return norm


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
	normnewfilex = norm(y, dx1vals, derivx(zvals1, x))
	normnewfiley = norm(y, dy1vals, derivy(zvals1, x))

	#norms grid2
	x1 = np.meshgrid(xyvals2, xyvals2)[0]
	y1 = np.meshgrid(xyvals2, xyvals2)[1]
	normnewfilex1 = norm(y1, dx2vals, derivx(zvals2, x1))
	normnewfiley1 = norm(y1, dy2vals, derivy(zvals2, x1))

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
	plt.figure()

	plt.scatter(x1, y1, c = normnewfilex1, marker= 'o',edgecolor='none')
	plt.colorbar()
	plt.show()

	plt.figure()

	plt.scatter(x1, y1, c = normnewfiley1, marker= 'o',edgecolor='none')
	plt.colorbar()
	plt.show()

	#e1 is the sum of all the norms for the values in the first file 


main()
