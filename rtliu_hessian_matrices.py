import yt
import numpy as np
import matplotlib.pyplot as plt
from math import log


def olr(a):  # outlier range for CB labelling
    q1 = np.percentile(a, 25)
    q3 = np.percentile(a, 75)
    iqr = q3-q1
    return([q1-1.5*iqr, q3+1.5*iqr])


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

    plt.xlim(0, 1)
    plt.ylim(-0.5, 0.5)
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


def quiverPlot(nx, ny, nc, cbname, fname, text=""):
    X, Y = nx, ny
    U = nc[0]
    V = nc[1]
    R = np.hypot(U, V)
    plt.clf()

    # plt.figure()
    # plt.title('Arrows scale with plot width, not view')
    # Q = plt.quiver(nx, ny, U, V, units='width')
    # qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
    #                    coordinates='figure')

    # plt.figure()
    # k = 2
    # Q = plt.quiver(X[::k, ::k], Y[::k, ::k], U[::k, ::k], V[::k, ::k], R,
    #                pivot='tip', units='xy', width=0.002, scale=1)
    # # qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
    # #                   coordinates='figure')
    # plt.scatter(X[::k, ::k], Y[::k, ::k], color='k', s=2)

    plt.quiver(X, Y, U, V, R, units='xy', angles='xy', pivot='tip',
               scale=1 / 0.1, minlength=0)
    # plt.scatter(X, Y, color='k', s=1)

    # plt.xlim(0, 1)
    # plt.ylim(-0.5,0.5)
    # plt.gca().set_position((.1, .3, .8, .6))

    plt.xlabel('x (cm)', fontsize=40)
    plt.ylabel('y (cm)', fontsize=40)
    plt.title(fname, fontsize=40)
    # cb = plt.colorbar()
    # cb.set_label(cbname, fontsize=40)
    plt.gcf().set_size_inches(60, 50)

    plt.savefig("{}.png".format(fname))


# Replaced with numpy array subtraction
def diff(a, b):  # Returns difference between datasets
    return([a[i] - b[i] for i in range(min(len(a), len(b)))])


def diffstats(ndiff):
    # Find relevant stats from ndiff
    # min = ndiff[0]
    # max = ndiff[0]
    # absmin = abs(ndiff[0])
    absmax = abs(float(ndiff[0][0]))
    RSS = 0
    for row in ndiff:
        for i in row:
            if abs(float(i)) > absmax:
                absmax = abs(float(i))

            RSS += float(i)**2

    return({"absmax": absmax, "absmaxr": round(absmax, 5),
           "RSS": RSS, "RSSr": round(RSS, 5)})


# Arguments are two datasets; specified parameter of interest;
# verbose parameter; color bar label; file name
def diffAnalysis(ds0, ds1, poi, poiv, cbname, fname):
    # ad0 = ds0.all_data()
    # ad1 = ds1.all_data()
    ad0 = ds0.covering_grid(level=0, left_edge=ds0.index.grids[0].LeftEdge,
                            dims=ds0.domain_dimensions)
    ad1 = ds1.covering_grid(level=0, left_edge=ds1.index.grids[0].LeftEdge,
                            dims=ds1.domain_dimensions)

    # Retrive relevant data from grid: x, y, poi
    nx = [np.array(ad0["x"]), np.array(ad1["x"])]
    ny = [np.array(ad0["y"]), np.array(ad1["y"])]

    npoi = [np.array(ad0[poi]), np.array(ad1[poi])]
    # ndiff = diff(npoi[0],npoi[1])
    ndiff = npoi[0]-npoi[1]

    stats = diffstats(ndiff)
    statsv = "Greatest Absolute Difference in {} = {}" \
        "\nSum of Squared Residuals (Measure of total discrepancy) = {}" \
        .format(poiv, stats["absmaxr"], stats["RSSr"])

    # plot(nx[0],ny[0],cbvely,nvely[0],"khv0") Uncomment asap
    # plot(nx[1],ny[1],cbvely,nvely[1],"khv1")
    plot(nx[0], ny[0], ndiff, cbname, fname, statsv)


# Arguments are a dataset; specified parameter of interest;
# verbose parameter; color bar label; file name
def visualize(ds, poi, poiv, cbname, fname):
    # ad = ds.all_data()
    ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                          dims=ds.domain_dimensions)

    # Retrive relevant data from grid: x, y, poi
    nx = np.array(ad["x"])
    ny = np.array(ad["y"])

    npoi = np.array(ad[poi])

    plot(nx, ny, npoi, cbname, fname)


# Arguments are two datasets; specified field of interest;
# verbose field; color bar label; file name
def fieldAnalysis(ds, foi, foiv, cbname, fname):
    ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                          dims=ds.domain_dimensions)

    x = np.array(ad["x"])
    y = np.array(ad["y"])

    plot(x, y, foi(ad), cbname, fname)


# Arguments are dataset; vector field;
# verbose field; color bar label; file name
def fieldAnalysisQuiver(ds, vf, foiv, cbname, fname):
    ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                          dims=ds.domain_dimensions)

    x = np.array(ad["x"])
    y = np.array(ad["y"])

    quiverPlot(x, y, vf, cbname, fname)


def vorticity(poi_0, poi_1):
    def f(ad):
        x = np.array(ad["x"])
        y = np.array(ad["y"])
        hx = float(x[1, 0]-x[0, 0])
        hy = float(y[0, 1]-y[0, 0])

        apoi_0 = np.array(ad[poi_0])
        apoi_1 = np.array(ad[poi_1])

        # dvelydx = np.zeros((len(x),len(x[0])))
        # dvelxdy = np.zeros((len(x),len(x[0])))
        dpoi_1dx = finiteDiffX(x, y, apoi_1)
        dpoi_0dy = finiteDiffY(x, y, apoi_0)
        vort = np.zeros((len(x), len(x[0])))

        # Finite differences used; with accuracy 2
        for i in range(len(x)):
            for j in range(len(x[i])):
                # Calculates vorticity
                vort[i, j] = dpoi_1dx[i, j] - dpoi_0dy[i, j]

                # Absolute Value
                # vort[i, j] = abs(dvelydx[i, j] - dvelxdy[i, j])

                # Log Scale
                # vort[i, j] = log(abs(dvelydx[i, j] - dvelxdy[i, j]))

        return(vort)

    return(f)


def jacobian(poi_0, poi_1):
    def f(ds):
        try:
            ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                                  dims=ds.domain_dimensions)
        except:
            ad = ds

        x = np.array(ad["x"])
        y = np.array(ad["y"])
        hx = float(x[1, 0]-x[0, 0])
        hy = float(y[0, 1]-y[0, 0])

        apoi_0 = np.array(ad[poi_0])
        apoi_1 = np.array(ad[poi_1])

        dpoi_0dx = finiteDiffX(x, y, apoi_0)
        dpoi_0dy = finiteDiffY(x, y, apoi_0)
        dpoi_1dx = finiteDiffX(x, y, apoi_1)
        dpoi_1dy = finiteDiffY(x, y, apoi_1)

        return([[dpoi_0dx, dpoi_0dy], [dpoi_1dx, dpoi_1dy]])

    return(f)


def hessian_old(ds, c):
    try:
        ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                              dims=ds.domain_dimensions)
    except:
        ad = ds

    x = np.array(ad["x"])
    y = np.array(ad["y"])

    ac = np.array(ad[c])

    dcdx2 = finiteDiffX2(x, y, ac)  # Second order derivative wrt x
    dcdy2 = finiteDiffY2(x, y, ac)  # Second order derivative wrt y
    dcdxy = finiteDiffY(x, y, finiteDiffX(x, y, ac))  # Mixed derivative wrt x y

    return([[dcdx2, dcdxy], [dcdxy, dcdy2]])


def hessian(ds, c):
    try:
        ad = ds.covering_grid(level=0, left_edge=ds.index.grids[0].LeftEdge,
                              dims=ds.domain_dimensions)
    except:
        ad = ds

    x = np.array(ad["x"])
    y = np.array(ad["y"])

    ac = np.array(ad[c])

    dcdx2 = finiteDiffX2(x, y, ac)  # Second order derivative wrt x
    dcdy2 = finiteDiffY2(x, y, ac)  # Second order derivative wrt y
    dcdxy = finiteDiffY(x, y, finiteDiffX(x, y, ac))  # Mixed derivative wrt x y

    vdy = finiteDiffX(x, y, ac)

    ah = [[[] for j in range(len(x[i]))] for i in range(len(x))]
    for i in range(len(x)):
        for j in range(len(x[i])):
            ah[i][j] = [[dcdx2[i][j], dcdxy[i][j]], [dcdxy[i][j], dcdy2[i][j]]]
            # if i == 0 and j < 100:
            #     print(ah[i][j])
                # print(dcdxy[i][j])
                # print(dcdx2[i][j])

    return(np.array(ah))


def dsDeterminant(a):
    ac = np.zeros((len(a), len(a[0])))
    for i in range(len(a)):
        for j in range(len(a[i])):
            m = a[i][j]
            det = m[0][0]*m[1][1] - m[0][1]*m[1][0]
            ac[i][j] = det
    return(ac)


def gradient(poi_0, poi_1, option="log"):
    def f(ad):
        j = jacobian(poi_0, poi_1)(ad)
        duxdx = j[0][0]
        duydy = j[1][1]
        gradMag = np.zeros((len(j[0][0]), len(j[0][0][0])))

        for i in range(len(gradMag)):
            for j in range(len(gradMag[i])):
                if option == "log":
                    gradMag[i, j] = log(1+(duxdx[i, j]**2 + duydy[i, j]**2)**(0.5))
                if option == "reg":
                    gradMag[i, j] = (duxdx[i, j]**2 + duydy[i, j]**2)**(0.5)

        return(gradMag)

    return(f)


def classify(poi):
    def f(ad):
        hy = hessian(ad, poi)
        return(dsDeterminant(hy))

    return(f)

def finiteDiffX(ax, ay, ac):
    x = np.array(ax)
    y = np.array(ay)
    c = np.array(ac)
    dc = np.zeros((len(x), len(x[0])))

    hx = float(x[1, 0]-x[0, 0])

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

    return(dc)


def finiteDiffY(ax, ay, ac):
    x = np.array(ax)
    y = np.array(ay)
    c = np.array(ac)
    dc = np.zeros((len(x), len(x[0])))

    hy = float(y[0, 1]-y[0, 0])

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


ds = [yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0000"),
      yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0001")]

# ds[0].print_stats()
# ds[1].print_stats()

# diffAnalysis(ds[0],ds[1],"velx","x-velocity","vel$_x$ (cm$\cdot$code length/code time)","KH_velx_analysis")
# diffAnalysis(ds[0],ds[1],"vely","y-velocity","vel$_y$ (cm$\cdot$code length/code time)","KH_vely_analysis")

# visualize(ds[0],"vely","y-velocity","vel$_y$ (cm$\cdot$code length/code time)","KH_vely_visualization_0")
# visualize(ds[1],"vely","y-velocity","vel$_y$ (cm$\cdot$code length/code time)","KH_vely_visualization_1")

# fieldAnalysis(ds[0],vorticity("velx", "vely"),"Vorticity","$\\vec{\\omega}$ (rad/second)","KH_vort_analysis_0")
# fieldAnalysis(ds[1],vorticity("velx", "vely"),"Vorticity","$\\vec{\\omega}$ (rad/second)","KH_vort_analysis_1")

# fieldAnalysis(ds[0],vorticity("velx", "vely"),"Vorticity","$\\vec{\\omega}$ (rad/second)","KH_vort_log_analysis_0")
# fieldAnalysis(ds[1],vorticity("velx", "vely"),"Vorticity","$\\vec{\\omega}$ (rad/second)","KH_vort_log_analysis_1")

# fieldAnalysis(ds[0],gradient("velx", "vely"),"Gradient","grad $\\vec{u}$","KH_grad_log_analysis2_0")
# fieldAnalysis(ds[1],gradient("velx", "vely"),"Gradient","grad $\\vec{u}$","KH_grad_log_analysis2_1")

# fieldAnalysisQuiver(ds[0],jacobian("velx", "vely")(ds[0])[0],"x-velocity Gradient","grad $\\vec{u}$","KH_grad_xquiv_analysis_")
# fieldAnalysisQuiver(ds[1],jacobian("velx", "vely")(ds[0])[1],"y-velocity Gradient","grad $\\vec{u}$","KH_grad_yquiv_analysis_")

fieldAnalysis(ds[0],gradient("magx", "magy", "reg"),"Gradient","grad $\\vec{u}$","KH_Bxygrad_analysis2_0_TEMP")
fieldAnalysis(ds[1],gradient("magx", "magy", "reg"),"Gradient","grad $\\vec{u}$","KH_Bxygrad_analysis2_1_TEMP")

# fieldAnalysis(ds[0],classify("vely"),"Determinant","det $\\vec{u}$","KH_Bxhessdet_analysis0_0")
# fieldAnalysis(ds[1],classify("vely"),"Determinant","det $\\vec{u}$","KH_Bxhessdet_analysis0_1")
