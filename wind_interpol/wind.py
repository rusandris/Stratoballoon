import ecmwf.data as ecdata
from pylab import *
from scipy.interpolate import interpn

#returns (u,v) wind velocity component from the dataset for requested (lat,lon,z)
def read_dataset(data, lat, lon, z, levels, lat_values, lon_values, reso=0.4, lat_min = -90,lat_max = 90,lon_min = -180,lon_max = 180):

    if lat not in lat_values or lon not in lon_values or z not in levels:
        return (None,None)
        
    uz = data.select(shortName = 'u',level=[z])
    vz = data.select(shortName = 'v',level=[z])

    grid =(int((lat_max - lat_min)/reso + 1), int((lon_max - lon_min)/reso))
    uz = reshape(uz.values(),grid)   
    vz = reshape(vz.values(),grid)
    
    i = int( np.round( (lat_max-lat)/reso , 1 ) )
    j = int( np.round( (lon-lon_min)/reso , 1 ) )
    #print(i,j)
    
    return (uz[i,j],vz[i,j])

#returns (u,v) wind velocity components (array of arrays of tuples) from the dataset for requested (lats,lons,zs)
def read_dataset_many(data, lats, lons, zs, levels, lat_values, lon_values, reso=0.4, lat_min = -90,lat_max = 90,lon_min = -180,lon_max = 180):
    
    return_set = []
    
    if len(lats) != len(lons):
        return "Dimension mismatch!"
    
    length = len(lats)
    
    for z in zs:
        set_for_z = []
        for i in range(length):
            
            lat = lats[i]
            lon = lons[i]

            if lat not in lat_values or lon not in lon_values or z not in levels:
                return (None,None)
            
                
            uz = data.select(shortName = 'u',level=[z])
            vz = data.select(shortName = 'v',level=[z])

            grid =(int((lat_max - lat_min)/reso + 1), int((lon_max - lon_min)/reso))
            uz = reshape(uz.values(),grid)   
            vz = reshape(vz.values(),grid)
            
            i = int( round( (lat_max-lat)/reso , 1 ) )
            j = int( round( (lon-lon_min)/reso , 1 ) )
            #print(i,j)
            
            set_for_z.append((uz[i,j],vz[i,j]))
        return_set.append(set_for_z)
    return return_set

#returns (u,v) wind velocity component for arbitrary coordinates (uses interpolation) 
#more interpol methods?
def get_wind_velocity(data,lat,lon,z, levels,lat_values,lon_values,interpol_method = "linear" ,reso=0.4, lat_min = -90,lat_max = 90,lon_min = -180,lon_max = 180):
    
    if interpol_method == "linear":
        return interpol_linear(data, lat, lon,z,levels,lat_values,lon_values)
    else:
        return "Interpolation method not supported"
        
def between(array,x):
    less = None
    greater = None
    i = 0
    while i < len(array) and x >= array[i]:
        i += 1
    if i != len(array):
        less = array[i-1]
        greater = array[i]
    elif x == array[-1]:
        less = array[i-2]
        greater = array[i-1]
        
    return (less,greater)

#linear interpolation on 3D grid
def interpol_linear(data, lat, lon,z,levels,lat_values,lon_values):

    neighbours_lat = between(lat_values[::-1],lat)
    neighbours_lon = between(lon_values,lon)
    neighbours_z = between(levels[::-1], z)
    #print(neighbours_lat)
    #print(neighbours_lon)
    #print(neighbours_z)

    if any([x is None for x in [*neighbours_lat,*neighbours_lon,*neighbours_z]]):
        return (None,None)
        
    grid_values = read_dataset_many(data,[neighbours_lat[0],neighbours_lat[0],neighbours_lat[1],neighbours_lat[1]],[neighbours_lon[0],neighbours_lon[1],neighbours_lon[0],neighbours_lon[1]],neighbours_z,levels,lat_values,lon_values)
    
    grid_values = array(grid_values)
    u_values = reshape(grid_values[:,:,0],(2,2,2))
    v_values = reshape(grid_values[:,:,1],(2,2,2))

    
    u = interpn((neighbours_z,neighbours_lat,neighbours_lon),u_values,[z,lat,lon])[0]
    v = interpn((neighbours_z,neighbours_lat,neighbours_lon),v_values,[z,lat,lon])[0]

    return (u,v)
    
#return (u,v) wind velocity components for all levels at (lat,lon)
def read_data_all_levels(data,lat,lon,levels,lat_values,lon_values):
    return array([read_dataset(data,lat,lon,z,levels,lat_values,lon_values) for z in levels])

#return (u,v) wind velocity components for all latitudes at (lon,z)
def read_data_loncircle(data,lon,z,levels,lat_values,lon_values):
    return array([read_dataset(data,lat,lon,z,levels,lat_values,lon_values) for lat in lat_values])
    
#return (u,v) wind velocity components for all longitudes at (lat,z)
def read_data_latcircle(data,lat,z,levels,lat_values,lon_values):
    return array([read_dataset(data,lat,lon,z,levels,lat_values,lon_values) for lon in lon_values])


            
