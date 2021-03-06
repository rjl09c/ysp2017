import yt
import numpy as np
import matplotlib.pyplot as plt
from math import log

def olr(a): #outlier range for CB labelling
	q1 = np.percentile(a,25)
	q3 = np.percentile(a,75)
	iqr = q3-q1
	return([q1-1.5*iqr,q3+1.5*iqr])

def arrMin(nc):
	try:
		return(nc.min())
	except:
		return(min(nc))

def arrMax(nc):
	try:
		return(nc.max())
	except:
		return(max(nc))

def plot(nx,ny,nc,cbname,fname,text=""):
	# print("{}.png loading".format(fname))
	print("{}: min = {} , max = {}".format(fname, arrMin(nc), arrMax(nc)))
	plt.clf()
	# plt.scatter(nx, ny, c=nc, marker='o', edgecolor='none', vmin=-0.01, vmax=0.01)

	# Uses extrema cutoffs for CB range
	plt.scatter(nx, ny, c=nc, marker='o', edgecolor='none', vmin=arrMin(nc), vmax=arrMax(nc))

	# Uses outlier cutoffs for CB range
	# plt.scatter(nx, ny, c=nc, marker='o', edgecolor='none', vmin=olr(nc)[0], vmax=olr(nc)[1])

	plt.xlim(0,1)
	plt.ylim(-0.5,0.5)
	plt.gca().set_position((.1, .3, .8, .6))

	plt.xlabel('x (cm)', fontsize=18)
	plt.ylabel('y (cm)', fontsize=18)
	plt.title(fname, fontsize=18)
	cb = plt.colorbar()
	cb.set_label(cbname, fontsize=18)
	plt.gcf().set_size_inches(12,10)

	plt.figtext(.02,.02,text)

	plt.savefig("{}.png".format(fname)) #Uncomment to save image
	# plt.show() #Uncomment to show image

def quiverPlot(nx,ny,nc,cbname,fname,text=""):
	X, Y = nx, ny
	U = nc[0]
	V = nc[1]

	plt.clf()

	# plt.figure()
	# plt.title('Arrows scale with plot width, not view')
	# Q = plt.quiver(nx, ny, U, V, units='width')
	# qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
	#                    coordinates='figure')

	plt.figure()
	plt.title(fname)
	k = 5
	Q = plt.quiver(X[::k, ::k], Y[::k, ::k], U[::k, ::k], V[::k, ::k], np.hypot(U, V),
	               pivot='tip')
	qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
	                   coordinates='figure')
	plt.scatter(X[::k, ::k], Y[::k, ::k], color='r', s=1)

	# plt.figure()
	# plt.title("pivot='tip'; scales with x view")
	# M = np.hypot(U, V)
	# Q = plt.quiver(X, Y, U, V, M, units='x', pivot='tip', width=0.001,
	#                scale=1 / 0.15)
	# plt.scatter(X, Y, color='k', s=1)	

	plt.savefig("{}.png".format(fname))

# Replaced with numpy array subtraction
def diff(a,b): #Returns difference between datasets
	return([a[i] - b[i] for i in range(min(len(a),len(b)))])

def diffstats(ndiff):
	# Find relevant stats from ndiff
	# min = ndiff[0]
	# max = ndiff[0]
	# absmin = abs(ndiff[0])
	absmax = abs(float(ndiff[0][0]))
	RSS = 0
	for row in ndiff:
		for i in row:
			# if i < velymin:
			# 	velymin = i
			# if i > velymax:
			# 	velymax = i
			# if abs(i) < velyabsmin:
			# 	velyabsmin = abs(i)
			if abs(float(i)) > absmax:
				absmax = abs(float(i))

			RSS += float(i)**2

	return({"absmax" : absmax, "absmaxr" : round(absmax,5), "RSS" : RSS, "RSSr" : round(RSS,5)})

# Arguments are two datasets; specified parameter of interest; verbose parameter; color bar label; file name
def diffAnalysis(ds0, ds1, poi, poiv, cbname, fname):
	# ad0 = ds0.all_data()
	# ad1 = ds1.all_data()
	ad0 = ds0.covering_grid(level=0, left_edge=ds0.index.grids[0].LeftEdge, dims=ds0.domain_dimensions)
	ad1 = ds1.covering_grid(level=0, left_edge=ds1.index.grids[0].LeftEdge, dims=ds1.domain_dimensions)
	
	# Retrive relevant data from grid: x, y, vely
	nx = [np.array(ad0["x"]), np.array(ad1["x"])]
	ny = [np.array(ad0["y"]), np.array(ad1["y"])]
	
	npoi = [np.array(ad0[poi]), np.array(ad1[poi])]
	# ndiff = diff(npoi[0],npoi[1])
	ndiff = npoi[0]-npoi[1]

	stats = diffstats(ndiff)
	statsv = "Greatest Absolute Difference in {} = {}" \
		"\nSum of Squared Residuals (Measure of total discrepancy) = {}" \
		.format(poiv,stats["absmaxr"],stats["RSSr"])

	# plot(nx[0],ny[0],cbvely,nvely[0],"khv0") Uncomment asap
	# plot(nx[1],ny[1],cbvely,nvely[1],"khv1")
	plot(nx[0],ny[0],ndiff,cbname,fname,statsv)

# Arguments are a dataset; specified parameter of interest; verbose parameter; color bar label; file name
def visualize(ds, poi, poiv, cbname, fname):
	# ad = ds.all_data()
	ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge, dims=ds.domain_dimensions)

	# Retrive relevant data from grid: x, y, vely
	nx = np.array(ad["x"])
	ny = np.array(ad["y"])
	
	npoi = np.array(ad[poi])

	plot(nx,ny,npoi,cbname,fname)

# Arguments are two datasets; specified field of interest; verbose field; color bar label; file name
def fieldAnalysis(ds, foi, foiv, cbname, fname):
	ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge, dims=ds.domain_dimensions)

	x = np.array(ad["x"])
	y = np.array(ad["y"])

	plot(x,y,foi(ad),cbname,fname)

def fieldAnalysisQuiver(ds, foi, foiv, cbname, fname):
	ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge, dims=ds.domain_dimensions)

	x = np.array(ad["x"])
	y = np.array(ad["y"])

	quiverPlot(x,y,foi(ad),cbname,fname)

def vorticity(ad):
	x = np.array(ad["x"])
	y = np.array(ad["y"])
	hx = float(x[1,0]-x[0,0])
	hy = float(y[0,1]-y[0,0])

	velx = np.array(ad["velx"])
	vely = np.array(ad["vely"])
	
	# dvelydx = np.zeros((len(x),len(x[0])))
	# dvelxdy = np.zeros((len(x),len(x[0])))
	dvelydx = finiteDiffX(x, y, vely)
	dvelxdy = finiteDiffY(x, y, velx)
	vort = np.zeros((len(x),len(x[0])))

	# Finite differences used; with accuracy 2
	for i in range(len(x)):
		for j in range(len(x[i])):
			# Calculates vorticity
			# vort[i,j] = dvelydx[i,j] - dvelxdy[i,j]

			# Absolute Value
			# vort[i,j] = abs(dvelydx[i,j] - dvelxdy[i,j])

			# Log Scale
			vort[i,j] = log(abs(dvelydx[i,j] - dvelxdy[i,j]))

			# Log Scale
			# vort[i,j] = log(abs(dvelydx[i,j]))

	return(vort)

def gradient(ad):
	x = np.array(ad["x"])
	y = np.array(ad["y"])
	hx = float(x[1,0]-x[0,0])
	hy = float(y[0,1]-y[0,0])

	velx = np.array(ad["velx"])
	vely = np.array(ad["vely"])
	
	# dvelydx = np.zeros((len(x),len(x[0])))
	# dvelxdy = np.zeros((len(x),len(x[0])))
	dvelxdx = finiteDiffX(x, y, velx)
	dvelydy = finiteDiffY(x, y, vely)
	gradMag = np.zeros((len(x),len(x[0])))

	# Finite differences used; with accuracy 2
	for i in range(len(x)):
		for j in range(len(x[i])):
			# Calculates vorticity
			# grad[i,j] = [dvelxdx[i,j], dvelydy[i,j]]

			# Magnitude
			# grad[i,j] = (dvelxdx[i,j]**2 + dvelydy[i,j]**2)**(0.5)

			# Log Scale
			gradMag[i,j] = log((dvelxdx[i,j]**2 + dvelydy[i,j]**2)**(0.5))

	# return(gradMag)
	return(dvelxdx, dvelydy)

def finiteDiffX(ax, ay, ac):
	x = np.array(ax)
	y = np.array(ay)
	c = np.array(ac)
	dc = np.zeros((len(x),len(x[0])))

	hx = float(x[1,0]-x[0,0])

	# Finite differences used; with accuracy 2
	for i in range(len(x)):
		for j in range(len(x[i])):
			if i > 0 and i < len(x)-1:
				dc[i,j] = (-0.5*c[i-1,j] + 0.5*c[i+1,j])/hx
			else:
				if i == 0:
					dc[i,j] = (-1.5*c[i,j] + 2*c[i+1,j] - 0.5*c[i+2,j])/hx
				if i == len(x)-1:
					dc[i,j] = (1.5*c[i-2,j] - 2*c[i-1,j] + 1.5*c[i,j])/hx

	return(dc)

def finiteDiffY(ax, ay, ac):
	x = np.array(ax)
	y = np.array(ay)
	c = np.array(ac)
	dc = np.zeros((len(x),len(x[0])))

	hy = float(y[0,1]-y[0,0])
	
	# Finite differences used; with accuracy 2
	for i in range(len(x)):
		for j in range(len(x[i])):
			if j > 0 and j < len(y)-1:
				dc[i,j] = (-0.5*c[i,j-1] + 0.5*c[i,j+1])/hy
			else:
				if j == 0:
					dc[i,j] = (-1.5*c[i,j] + 2*c[i,j+1] - 0.5*c[i,j+2])/hy
				if j == len(y)-1:
					dc[i,j] = (1.5*c[i,j-2] - 2*c[i,j-1] + 1.5*c[i,j])/hy

	return(dc)

ds = [yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0000"), \
	yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0001")]

# ds[0].print_stats()
# ds[1].print_stats()

# diffAnalysis(ds[0],ds[1],"velx","x-velocity","vel$_x$ (cm$\cdot$code length/code time)","KH_velx_analysis")
# diffAnalysis(ds[0],ds[1],"vely","y-velocity","vel$_y$ (cm$\cdot$code length/code time)","KH_vely_analysis")

# visualize(ds[0],"vely","y-velocity","vel$_y$ (cm$\cdot$code length/code time)","KH_vely_visualization_0")
# visualize(ds[1],"vely","y-velocity","vel$_y$ (cm$\cdot$code length/code time)","KH_vely_visualization_1")

# fieldAnalysis(ds[0],vorticity,"Vorticity","$\\vec{\\omega}$ (rad/second)","KH_vort_log_analysis_0")
# fieldAnalysis(ds[1],vorticity,"Vorticity","$\\vec{\\omega}$ (rad/second)","KH_vort_log_analysis_1")

# fieldAnalysis(ds[0],gradient,"Gradient","grad $\\vec{u}$","KH_grad_log_analysis_0")
# fieldAnalysis(ds[1],gradient,"Gradient","grad $\\vec{u}$","KH_grad_log_analysis_1")

fieldAnalysisQuiver(ds[0],gradient,"Gradient","grad $\\vec{u}$","KH_grad_quiv_analysis_0")
fieldAnalysisQuiver(ds[1],gradient,"Gradient","grad $\\vec{u}$","KH_grad_quiv_analysis_1")