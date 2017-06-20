import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab
from yt.analysis_modules.halo_finding.api import HaloFinder


ds = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0000")
ds1 = yt.load("kh_mhd_Ma=0.803333333333At=0.0hdf5_chk_0001")

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

plt.scatter(xcoords, ycoords,c=vely-vely1, marker= 'o',edgecolor='none')
cb = plt.colorbar()
cb.set_label('Difference in y Velocity')
s=yt.SlicePlot(ds, 'z', ('gas','velocity_y'), center='m')
s.save()

plt.scatter(xcoords1, ycoords1,c=vely-vely1, marker= 'o',edgecolor='none')
cb = plt.colorbar()
cb.set_label('Difference in y Velocity')
r=yt.SlicePlot(ds1, 'z', ('gas','velocity_y'), center='m')
r.save()

sumofy=0
for i in range(len(ycoords)):
	sumofy=sumofy+ ((vely[i]-vely1[i])**2)

sumofy=sumofy**(.5)

print(sumofy)
