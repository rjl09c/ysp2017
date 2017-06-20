import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder
import h5py

def vorticity(velx):
	vorticity=np.zeros((320,320))
	x2=0
	x3=0
	y1=0
	y2=0
	y3=0
	for i in range(len(velx)):
		for x in range(len(velx)):
			if 0<i<len(velx)-1 and 0<x<len(velx)-1:
				x1=velx[i-1][x]
				x3=velx[i+1][x]
				y1=velx[i][x-1]
				y3=velx[i][x+1]
				vorticity[i,x]=(((((-1/2)*x1)+((1/2)*x3))) - (((-1/2)*y1)+((1/2)*y3)))
			elif i==0 and 0<x<len(velx)-1:
				x1=velx[i][x]
				x2=velx[i+1][x]
				x3=velx[i+2][x]
				y1=velx[i][x-1]
				y2=velx[i][x]
				y3=velx[i][x+1]
				vorticity[i,x]=((((-3/2)*x1)+(2*x2)+((-1/2)*x3)) - (((-3/2)*y1)+(2*y2)+((-1/2)*y3)))
			elif x==0 and 0<i<len(velx)-1:
				x1=velx[i-1][x]
				x3=velx[i+1][x]
				y1=velx[i][x]
				y2=velx[i][x+1]
				y3=velx[i][x+2]
				vorticity[i,x]=(((((-1/2)*x1)+((1/2)*x3))) - (((-3/2)*y1)+(2*y2)+((-1/2)*y3)))
			elif i==0 and x==0:
				x1=velx[i][x]
				x2=velx[i+1][x]
				x3=velx[i+2][x]
				y1=velx[i][x]
				y2=velx[i][x+1]
				y3=velx[i][x+2]
				vorticity[i,x]=((((-3/2)*x1)+(2*x2)+((-1/2)*x3)) - (((-3/2)*y1)+(2*y2)+((-1/2)*y3)))
			elif i==len(velx)-1 and 0<x<len(velx)-1:
				x1=velx[i][x]
				x2=velx[i-1][x]
				x3=velx[i-2][x]
				y1=velx[i][x-1]
				y3=velx[i][x+1]
				vorticity[i,x]=((((-3/2)*x1)+(2*x2)+((-1/2)*x3) - (((-1/2)*y1)+((1/2)*y3))))
			elif x==len(velx)-1 and 0<i<len(velx)-1:
				x1=velx[i-1][x]
				x3=velx[i+1][x]
				y1=velx[i][x]
				y2=velx[i][x-1]
				y3=velx[i][x-2]
				vorticity[i,x]=(((((-1/2)*x1)+((1/2)*x3))) - (((-3/2)*y1)+(2*y2)+((-1/2)*y3)))
			elif x==len(velx)-1 and i==len(velx)-1:
				x1=velx[i][x]
				x2=velx[i-1][x]
				x3=velx[i-2][x]
				y1=velx[i][x]
				y2=velx[i][x-1]
				y3=velx[i][x-2]
				vorticity[i,x]=((((-3/2)*x1)+(2*x2)+((-1/2)*x3) - (((-3/2)*y1)+(2*y2)+((-1/2)*y3))))
		#x=0
	return vorticity


def finalvorticity(xcoords,velx, vely):
	vorticityx=vorticity(velx)
	vorticityy=vorticity(vely)
	vorticityfinal=np.zeros((320,320))
	distance=xcoords[1][0]-xcoords[0][0]
	for i in range(len(vorticityx)):
		for x in range(len(vorticityx)):
			vorticityfinal[i][x]=(vorticityx[i][x]-vorticityy[i][x])/distance
		#x=0
	return vorticityfinal

ds = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0000")
ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge, 
									dims=ds.domain_dimensions)

xcoords=ad["x"]
ycoords=ad["y"]
velx=ad["velx"]
vely=ad["vely"]

'''#writes vorticity data to a new file
f = h5py.File("vorticity.h5", "w")
f.create_dataset("/vorticity", data=vorticity(velx,vely))
f.close()
f = h5py.File("vorticity.h5", "r")
print(f["vorticity"].value)'''

ds1 = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0001")
dd = ds1.covering_grid(level=0, left_edge=ds1.index.grids[0].LeftEdge, dims=ds1.domain_dimensions)

xcoords1=dd["x"]
ycoords1=dd["y"]
velx1=dd["velx"]
vely1=dd["vely"]

'''
#writes vorticity data for second file to a new file
#f1 = h5py.File("vorticity1.h5", "w")
#f1.create_dataset("/vorticity1", data=vorticity(velx1,vely1))
#f1.close()
#f1 = h5py.File("vorticity1.h5", "r")
#print(f["vorticity1"].value)
'''


plt.scatter(xcoords, ycoords,c=finalvorticity(xcoords,velx,vely), marker= 'o',edgecolor='none')
cb = plt.colorbar()
cb.set_label('Vorticity')
plt.show()
#s=yt.SlicePlot(ds, 'z', ('gas','velocity_y') , center ='m')
#s.save()

plt.scatter(xcoords1, ycoords1,c=finalvorticity(xcoords1,velx1,vely1), marker= 'o',edgecolor='none')
cb = plt.colorbar()
cb.set_label('Vorticity')
plt.show()
#r=yt.SlicePlot(ds1, 'z', ('gas','vorticity_y'), center ='m')
#r.save()

	
