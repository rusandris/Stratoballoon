import ecmwf.data as ecdata
from pylab import *
from scipy.interpolate import interpn
from wind import *
import time


data = ecdata.read("20220725000000-3h-oper-fc.grib2")
levels = [1000, 925, 850, 700, 500, 300, 250, 200, 50]

reso=0.4
lat_min = -90
lat_max = 90
lon_min = -180
lon_max = 180

# latitude and longitude from dataset
lat_values = np.round(data.latitudes()[0][0::int((lon_max - lon_min)/reso)],1)
lon_values = np.round(data.longitudes()[0][0:int((lon_max - lon_min)/reso)],1)

#print(read_dataset(data,25.2,36.4,700,levels))

################################### Z interpolation #######################################
ztest_data = read_data_all_levels(data,24.8,36.4,levels,lat_values,lon_values)
#print(test_data)
us = ztest_data[:,0]
vs = ztest_data[:,1]

zs = arange(51,1000,20)

zinterp = array([get_wind_velocity(data,24.8,36.4,z,levels,lat_values,lon_values) for z in zs])
zinterp_u = zinterp[:,0]
zinterp_v = zinterp[:,1]

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Z interpolation test')

ax1.plot(levels,us,'.-',c='k')
ax1.scatter(zs,zinterp_u,c='r',marker='x')
ax1.set(xlabel="z", ylabel="u (m/s)")



ax2.plot(levels,vs,'.-',c='k')
ax2.scatter(zs,zinterp_v,c='r',marker='x')
ax2.set(xlabel="z", ylabel="v (m/s)")
show()


################################### angle interpolation #######################################
lattest_data = read_data_loncircle(data,36.4,1000,levels,lat_values,lon_values)
us = lattest_data[:,0]
vs = lattest_data[:,1]


lats = arange(-90,91,0.2)

latinterp = array([get_wind_velocity(data,lat,36.4,1000,levels,lat_values,lon_values) for lat in lats])
latinterp_u = latinterp[:,0]
latinterp_v = latinterp[:,1]




fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('lat interpolation test')

ax1.plot(lat_values,us,'.-',c='k')
ax1.scatter(lats,latinterp_u,c='r',marker='x')
ax1.set(xlabel="latitude", ylabel="u (m/s)")


ax2.plot(lat_values,vs,'.-',c='k')
ax2.scatter(lats,latinterp_v,c='r',marker='x')
ax2.set(xlabel="latitude", ylabel="v (m/s)")
show()


################################### angle interpolation #######################################
lontest_data = read_data_latcircle(data,24.8,1000,levels,lat_values,lon_values)
us = lontest_data[:,0]
vs = lontest_data[:,1]


lons = arange(-180,179.6,0.2)

loninterp = array([get_wind_velocity(data,24.8,lon,1000,levels,lat_values,lon_values) for lon in lons])
loninterp_u = loninterp[:,0]
loninterp_v = loninterp[:,1]


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('lon interpolation test')

ax1.plot(lon_values,us,'.-',c='k')
ax1.scatter(lons,loninterp_u,c='r',marker='x')
ax1.set(xlabel="longitude", ylabel="u (m/s)")


ax2.plot(lon_values,vs,'.-',c='k')
ax2.scatter(lons,loninterp_v,c='r',marker='x')
ax2.set(xlabel="longitude", ylabel="v (m/s)")
show()




