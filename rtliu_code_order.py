import numpy as np
import matplotlib.pyplot as plt
from math import log


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


def plot(nx, ny, nc, cbname, fname, text=""):
    # print("{}.png loading".format(fname))
    print("{}: min = {} , max = {}".format(fname, arrMin(nc), arrMax(nc)))
    plt.clf()
    # plt.scatter(nx, ny, c=nc, marker='o', vmin=-0.01, vmax=0.01)

    # Uses extrema cutoffs for CB range
    plt.scatter(nx, ny, c=nc, marker='o', vmin=arrMin(nc), vmax=arrMax(nc))

    # Uses outlier cutoffs for CB range
    # plt.scatter(nx, ny, c=nc, marker='o', vmin=olr(nc)[0], vmax=olr(nc)[1])

    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.gca().set_position((.1, .3, .8, .6))

    plt.xlabel('x (cm)', fontsize=18)
    plt.ylabel('y (cm)', fontsize=18)
    plt.title(fname, fontsize=18)
    cb = plt.colorbar()
    cb.set_label(cbname, fontsize=18)
    plt.gcf().set_size_inches(12, 10)

    plt.figtext(.02, .02, text)

    plt.savefig("{}.png".format(fname))  # Uncomment to save image
    # plt.show() #Uncomment to show image


def finiteDiffX(ax, ay, ac):
    x = np.array(ax)
    y = np.array(ay)
    c = np.array(ac)
    dc = np.zeros((len(x), len(x[0])))
    hx = float(x[1, 0] - x[0, 0])
    # Finite differences used; with accuracy 2
    for i in range(len(x)):
        for j in range(len(x[i])):
            if i > 0 and i < len(x)-1:
                dc[i, j] = (-0.5*c[i-1, j] + 0.5*c[i+1, j])/hx
            else:
                if i == 0:
                    dc[i, j] = (-1.5*c[i, j] + 2*c[i+1, j] - 0.5*c[i+2, j])/hx
                if i == len(x)-1:
                    dc[i, j] = (0.5*c[i-2, j] - 2*c[i-1, j] + 1.5*c[i, j])/hx
                # dc[i, j] = 0
                # if i == 0:
                #     dc[i, j] = (-11/6*c[i, j] + 3*c[i+1, j] - 1.5*c[i+2, j] + 1/3*c[i+3, j])/hx
                # if i == len(x)-1:
                #     dc[i, j] = (-1/3*c[i-3, j] + 1.5*c[i-2, j] - 3*c[i-1, j] + 11/6*c[i, j])/hx
    return(dc)


def finiteDiffY(ax, ay, ac):
    x = np.array(ax)
    y = np.array(ay)
    c = np.array(ac)
    dc = np.zeros((len(x), len(x[0])))

    hy = float(y[0, 1] - y[0, 0])

    # Finite differences used; with accuracy 2
    for i in range(len(x)):
        for j in range(len(x[i])):
            if j > 0 and j < len(y)-1:
                dc[i, j] = (-0.5*c[i, j-1] + 0.5*c[i, j+1])/hy
            else:
                if j == 0:
                    dc[i, j] = (-1.5*c[i, j] + 2*c[i, j+1] - 0.5*c[i, j+2])/hy
                if j == len(y)-1:
                    dc[i, j] = (0.5*c[i, j-2] - 2*c[i, j-1] + 1.5*c[i, j])/hy
                # dc[i, j] = 0
                # if j == 0:
                #     dc[i, j] = (-11/6*c[i, j] + 3*c[i, j+1] - 1.5*c[i, j+2] + 1/3*c[i, j+3])/hy
                # if j == len(x)-1:
                #     dc[i, j] = (-1/3*c[i, j-3] + 1.5*c[i, j-2] - 3*c[i, j-1] + 11/6*c[i, j])/hy

    return(dc)


# Second derivative with respect to x
def finiteDiffX2(ax, ay, ac):
    x = np.array(ax)
    y = np.array(ay)
    c = np.array(ac)
    dc = np.zeros((len(x), len(x[0])))

    hx = float(x[1, 0]-x[0, 0])

    # Finite differences used; with accuracy 2
    for i in range(len(x)):
        for j in range(len(x[i])):
            if i > 0 and i < len(x)-1:
                dc[i, j] = (c[i-1, j] - 2*c[i, j] + c[i+1, j])/hx
            else:
                if i == 0:
                    dc[i, j] = (2*c[i, j] - 5*c[i+1, j] +
                                4*c[i+2, j] - c[i+3, j])/hx
                if i == len(x)-1:
                    dc[i, j] = (-1*c[i-3, j] + 4*c[i-2, j] -
                                5*c[i-1, j] + 2*c[i, j])/hx
                # dc[i, j] = 0

    return(dc)


# Second derivative with respect to y
def finiteDiffY2(ax, ay, ac):
    x = np.array(ax)
    y = np.array(ay)
    c = np.array(ac)
    dc = np.zeros((len(x), len(x[0])))

    hy = float(y[0, 1]-y[0, 0])

    # Finite differences used; with accuracy 2
    for i in range(len(x)):
        for j in range(len(x[i])):
            if j > 0 and j < len(y)-1:
                dc[i, j] = (c[i, j-1] - 2*c[i, j] + c[i, j+1])/hy
            else:
                if j == 0:
                    dc[i, j] = (2*c[i, j] - 5*c[i, j+1] +
                                4*c[i, j+2] - c[i, j+3])/hy
                if j == len(y)-1:
                    dc[i, j] = (-1*c[i, j-3] + 4*c[i, j-2] -
                                5*c[i, j-1] + 2*c[i, j])/hy
                # dc[i, j] = 0

    return(dc)


# Arguments are a dataset; specified field of interest;
# verbose field; color bar label; file name
def fieldAnalysis(ds, foi, foiv, cbname, fname):
    try:
        ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                              dims=ds.domain_dimensions)
    except:
        ad = ds

    x = np.array(ad["x"])
    y = np.array(ad["y"])

    plot(x, y, foi(ad), cbname, fname)


# Arguments are a dataset; specified field of interest;
# verbose field; color bar label; file name
def derivAnalysis(ds, foi, foiv, cbname, fname):
    try:
        ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                              dims=ds.domain_dimensions)
    except:
        ad = ds

    x = np.array(ad["x"])
    y = np.array(ad["y"])
    z = np.array(ad["z"])

    plot(x, y, foi(x, y, z), cbname, fname)


# Arguments are a dataset; specified parameter of interest;
# verbose parameter; color bar label; file name
def visualize(ds, poi, poiv, cbname, fname):
    try:
        # ad = ds.all_data()
        ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                              dims=ds.domain_dimensions)
    except:
        ad = ds

    # Retrive relevant data from grid: x, y, poi
    nx = np.array(ad["x"])
    ny = np.array(ad["y"])

    npoi = np.array(ad[poi])

    plot(nx, ny, npoi, cbname, fname)


# Handles csv file for compatability
def handle(fnameXY, fnameZ):
    fXY = open(fnameXY)
    fZ = open(fnameZ)
    grid = list(map(float, fXY.readline().split(",")))
    gridZ = [list(map(float, i.replace("\n", "").split(","))) for i in fZ.readlines()]
    fXY.close()
    fZ.close()
    x = [[i for j in range(len(grid))] for i in grid]
    y = [[j for j in grid] for i in range(len(grid))]
    return({"x": x, "y": y, "z": gridZ})


def normstats(ndiff):
    # Find relevant stats from ndiff
    # min = ndiff[0]
    # max = ndiff[0]
    # absmin = abs(ndiff[0])
    absmax = abs(float(ndiff[0][0]))
    err = 0
    for row in ndiff:
        for i in row:
            if abs(float(i)) > absmax:
                absmax = abs(float(i))

            err += abs(float(i))

    return({"absmax": absmax, "absmaxr": round(absmax, 10),
           "err": err, "errr": round(err, 10), "avgerr": round(err/(len(ndiff)**2), 10)})


# Arguments are two datasets; specified parameter of interest;
# verbose parameter; color bar label; file name
def normAnalysis(ds, ds_true, f, poiv, cbname, fname):
    try:
        # ad0 = ds0.all_data()
        # ad1 = ds1.all_data()
        ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                                dims=ds.domain_dimensions)
        ad_true = ds_true.covering_grid(level=0, left_edge=ds_true.index.grids[0].LeftEdge,
                                dims=ds_true.domain_dimensions)
    except:
        ad = ds
        ad_true = ds_true

    # Retrive relevant data from grid: x, y, poi
    nx = [np.array(ad["x"]), np.array(ad_true["x"])]
    ny = [np.array(ad["y"]), np.array(ad_true["y"])]

    npoi = [np.array(f(ad["x"], ad["y"], ad["z"])), np.array(ad_true["z"])]
    ndiff = npoi[0] - npoi[1]

    stats = normstats(ndiff)
    statsv = "Greatest Absolute Difference in {} = {}" \
        "\Average of Differences (Measure of total difference) = {}" \
        .format(poiv, stats["absmaxr"], stats["avgerr"])

    plot(nx[0], ny[0], ndiff, cbname, fname, statsv)

    return(stats["avgerr"])


ds = [handle("xy1.csv", "grid1.csv"), handle("xy2.csv", "grid2.csv")]
dsY = [handle("xy1.csv", "Dx1.csv"), handle("xy2.csv", "Dx2.csv")]
dsX = [handle("xy1.csv", "Dy1.csv"), handle("xy2.csv", "Dy2.csv")]

visualize(ds[0],"z","z values","$z$ (cm)","data1_z_visualization_a")
# visualize(ds[1],"z","z values","$z$ (cm)","data2_z_visualization_a")

derivAnalysis(ds[0], finiteDiffX, "$\\frac{\\partial }{\\partial x}$ z", "$\\frac{\\partial }{\\partial x}$ z","data1_ddx1_analysis_a")
# derivAnalysis(ds[0], finiteDiffY, "$\\frac{\\partial }{\\partial y}$ z", "$\\frac{\\partial }{\\partial y}$ z","data1_ddy1_analysis_a")

# derivAnalysis(ds[0], finiteDiffX, "$\\frac{\\partial^2 }{\\partial x^2}$ z", "$\\frac{\\partial^2 }{\\partial x^2}$ z","data1_ddx1_analysis_a")
# derivAnalysis(ds[0], finiteDiffY, "$\\frac{\\partial^2 }{\\partial y^2}$ z", "$\\frac{\\partial^2 }{\\partial y^2}$ z","data1_ddy1_analysis_a")

visualize(dsX[0],"z","dx values","dx","data1_dx_visualization_a")
# visualize(dsY[0],"z","dy values","dx","data1_dy_visualization_a")



edx_1 = normAnalysis(ds[0],dsX[0],finiteDiffX,"z diff","$\\frac{\\partial }{\\partial x}$ z","data1_zdiffdx_analysis_b")
edy_1 = normAnalysis(ds[0],dsY[0],finiteDiffY,"z diff","$\\frac{\\partial }{\\partial y}$ z","data1_zdiffdy_analysis_b")

edx_2 = normAnalysis(ds[1],dsX[1],finiteDiffX,"z diff","$\\frac{\\partial }{\\partial x}$ z","data2_zdiffdx_analysis_b")
edy_2 = normAnalysis(ds[1],dsY[1],finiteDiffY,"z diff","$\\frac{\\partial }{\\partial y}$ z","data2_zdiffdy_analysis_b")

h1 = (-1.9596) - (-2)
h2 = (-1.9799) - (-2)
print("Slope for dx: {}".format(log(edx_1/edx_2)/log(h1/h2)))
print("Slope for dy: {}".format(log(edy_1/edy_2)/log(h1/h2)))