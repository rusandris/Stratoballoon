from pylab import *
from wind_velocity import *

data = ecdata.read("20220725000000-3h-oper-fc.grib2")
levels = [1000, 925, 850, 700, 500, 300, 250, 200, 50]
lat_lims = (40,50)
lon_lims = (20,30)

dataf,lats,lons = filter_to_region(data,lat_lims,lon_lims,levels)


########################### height interpolation test #####################
zs = arange(50,1000,20)

us = [list(dataf.values())[i][0][0][0] for i in range(len(levels))]
vs = [list(dataf.values())[i][1][0][0] for i in range(len(levels))]

zinterp = array([interpol_linear(dataf,lats[0],lons[0],z,levels,lats,lons) for z in zs])
zinterp_u = zinterp[:,0]
zinterp_v = zinterp[:,1]

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Z interpolation test')

ax1.plot(levels[::-1],us,'.-',c='k')
ax1.scatter(zs,zinterp_u,c='r',marker='x')
ax1.set(xlabel="z (mbar)", ylabel="u (m/s)")

ax2.plot(levels[::-1],vs,'.-',c='k')
ax2.scatter(zs,zinterp_v,c='r',marker='x')
ax2.set(xlabel="z (mbar)", ylabel="v (m/s)")
show()


######################## longitude interpolation test #######################

longitude = arange(20,30,0.1)

us = list(dataf.values())[-1][0][0] 
vs = list(dataf.values())[-1][1][0]

loninterp = array([interpol_linear(dataf,lats[0],lon,1000,levels,lats,lons) for lon in longitude])
loninterp_u = loninterp[:,0]
loninterp_v = loninterp[:,1]


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('lon interpolation test')

ax1.plot(lons,us,'.-',c='k')
ax1.scatter(longitude,loninterp_u,c='r',marker='x')
ax1.set(xlabel="longitude", ylabel="u (m/s)")

ax2.plot(lons,vs,'.-',c='k')
ax2.scatter(longitude,loninterp_v,c='r',marker='x')
ax2.set(xlabel="longitude", ylabel="v (m/s)")

show()



######################## latitude interpolation test #######################

latitude = arange(40,50,0.1)

us = [list(dataf.values())[-1][0][i][0] for i in range(len(lats))]
vs = [list(dataf.values())[-1][1][i][0] for i in range(len(lats))]


latinterp = array([interpol_linear(dataf,lat,lons[0],1000,levels,lats,lons) for lat in latitude])
latinterp_u = latinterp[:,0]
latinterp_v = latinterp[:,1]


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('lat interpolation test')

ax1.plot(lats,us,'.-',c='k')
ax1.scatter(latitude,latinterp_u,c='r',marker='x')
ax1.set(xlabel="latitude", ylabel="u (m/s)")

ax2.plot(lats,vs,'.-',c='k')
ax2.scatter(latitude,latinterp_v,c='r',marker='x')
ax2.set(xlabel="latitude", ylabel="v (m/s)")

show()

