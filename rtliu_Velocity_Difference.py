import yt
import numpy as np
import matplotlib.pyplot as plt

def olr(a): #outlier range for CB labelling
	q1 = np.percentile(a,25)
	q3 = np.percentile(a,75)
	iqr = q3-q1
	return([q1-1.5*iqr,q3+1.5*iqr])

def plot(nx,ny,nc,cbname,fname,text=""):
	print("{}: min = {} , max = {}".format(fname, min(nc), max(nc)))
	plt.clf()
	# plt.scatter(nx, ny, c=nc, marker='o', edgecolor='none', vmin=-0.01, vmax=0.01)
	plt.scatter(nx, ny, c=nc, marker='o', edgecolor='none', vmin=olr(nc)[0], vmax=olr(nc)[1])
	plt.xlim(0,1)
	plt.ylim(-0.5,0.5)
	plt.gca().set_position((.1, .3, .8, .6))

	plt.xlabel('x (cm)')
	plt.ylabel('y (cm)')
	plt.title(fname)
	cb = plt.colorbar()
	cb.set_label(cbname)
	plt.gcf().set_size_inches(12,10)

	plt.figtext(.02,.02,text)

	plt.savefig("{}.png".format(fname)) #Uncomment to save image
	# plt.show() #Uncomment to show image

def diff(a,b): #Returns difference between datasets
	return([a[i] - b[i] for i in range(min(len(a),len(b)))])

def diffstats(ndiff):
	# Find relevant stats from ndiff
	# min = ndiff[0]
	# max = ndiff[0]
	# absmin = abs(ndiff[0])
	absmax = abs(ndiff[0])
	RSS = 0
	for i in ndiff:
		# if i < velymin:
		# 	velymin = i
		# if i > velymax:
		# 	velymax = i
		# if abs(i) < velyabsmin:
		# 	velyabsmin = abs(i)
		if abs(i) > absmax:
			absmax = abs(i)

		RSS += i**2

	return({"absmax" : absmax, "absmaxr" : round(absmax,5), "RSS" : RSS, "RSSr" : round(RSS,5)})


# Arguments are two datasets; specified parameter of interest; verbose parameter; color bar label; file name
def diffAnalysis(ds0, ds1, poi, poiv, cbname, fname):
	ad0 = ds0.all_data()
	ad1 = ds1.all_data()
	
	# Retrive relevant data from grid: x, y, vely
	nx = [np.array(ad0["x"]), np.array(ad1["x"])]
	ny = [np.array(ad0["y"]), np.array(ad1["y"])]
	
	npoi = [np.array(ad0[poi]), np.array(ad1[poi])]
	ndiff = diff(npoi[0],npoi[1])

	stats = diffstats(ndiff)
	statsv = "Greatest Absolute Difference in {} = {}\n\
		Sum of Squared Residuals (Measure of total discrepancy) = {}"\
		.format(poiv,stats["absmaxr"],stats["RSSr"])

	# plot(nx[0],ny[0],cbvely,nvely[0],"khv0") Uncomment asap
	# plot(nx[1],ny[1],cbvely,nvely[1],"khv1")
	plot(nx[0],ny[0],ndiff,cbname,fname,statsv)

ds = [yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0000"), \
	yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0001")]

ds[0].print_stats()
ds[1].print_stats()
# print(ds0.derived_field_list)
# print(ds1.derived_field_list)

# print(ds0.field_list)
# # print(ds1.field_list)
# # print(ds0.index.grid_left_edge)
# g0 = ds0.index.grids[1]
# g1 = ds1.index.grids[1]
# print(g0["velz"][:,:,0])
# print(g1["velz"][:,:,0])

diffAnalysis(ds[0],ds[1],"velx","x-velocity","vel$_x$ (cm$\cdot$code length/code time","KH_velx_analysis")
diffAnalysis(ds[0],ds[1],"vely","y-velocity","vel$_y$ (cm$\cdot$code length/code time","KH_vely_analysis")