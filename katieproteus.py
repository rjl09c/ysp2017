import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder
import os

#current files are from 20170706_1407
#previous files are from 20170705_0400

#finds the file
def find(name, path):
	for root, dirs, files in os.walk(path):
		if name in files:
			return os.path.join(root, name)


#plotting function
def plot(x, y, diff, fileEndName, colorname, plotTitle):
	plt.figure()
	plt.scatter(x, y,c=diff, marker= 'o',edgecolor='none')
	plt.title(plotTitle)
	cb = plt.colorbar()
	cb.set_label(colorname)
	plt.savefig(fileEndName+".png")


#calculating l1 norm 
def  norm(ycoords, vely, vely1):
	sumofy = 0
	for i in range(len(ycoords)):
		sumofy=sumofy+ ((vely[i]-vely1[i])**2)

	sumofy=sumofy**(.5)
	return sumofy

#get 1st and last files and retrieve info
def get_files(file1current, file1previous, num, title):
	ds = yt.load(file1current)
	ds1 = yt.load(file1previous)
	
	ad=ds.all_data()
	xcoords=ad["x"]
	ycoords=ad["y"]
	velx=ad["velx"]
	vely=ad["vely"]

	dd=ds1.all_data()
	xcoords1=dd["x"]
	ycoords1=dd["y"]
	velx1=dd["velx"]
	vely1=dd["vely"]

	plot(xcoords,ycoords,velx - velx1, "chk_00" + num + "x" , "Velocity X Diff", "Velocity X Analysis" + title)
	plot(xcoords1,ycoords1,vely - vely1, "chk_00" + num, "Velocity Y Diff", "Velocity Y Analysis" + title)

	print(norm(ycoords, velx, velx1))
	print(norm(ycoords, vely, vely1))


def main():
	'''when using server
	get_files(find("sedov_2d_4lev_hdf5_chk_0000", "~/yspKatherine/ysp/ProteusTest/results/20170706_1407/comparison/suite/sedov/sedov_2d_4lev"), 
				find("sedov_2d_4lev_hdf5_chk_0000", "~/yspKatherine/ysp/ProteusTest/results/20170706_1007/comparison/suite/sedov/sedov_2d_4lev"),
						 "00", "for 1st Checkpoint")

	get_files(find("sedov_2d_4lev_hdf5_chk_0056", "~/yspKatherine/ysp/ProteusTest/results/20170705_0400/comparison/suite/sedov/sedov_2d_4lev"), 
				find("sedov_2d_4lev_hdf5_chk_0056", "~/yspKatherine/ysp/ProteusTest/results/20170705_0400/comparison/suite/sedov/sedov_2d_4lev"),
						 "56", "for 2nd Checkpoint")
	'''

	get_files("sedov_2d_4lev_hdf5_chk_0000", "sedov1_2d_4lev_hdf5_chk_0000","00", " for 0000 Checkpoint")

	get_files("sedov_2d_4lev_hdf5_chk_0056", "sedov1_2d_4lev_hdf5_chk_0056","56", " for 0056 Checkpoint")


main()

