import ecmwf.data as ecdata
from pylab import *
from scipy.interpolate import interpn


#filters dataset to region -> returns dictionary structure z : (u_matrix,v_matrix)  
def filter_to_region(data,lat_lims,lon_lims,levels,reso=0.4, lat_min = -90,lat_max = 90,lon_min = -180,lon_max = 180):
    lat_low,lat_high = lat_lims
    lon_low,lon_high = lon_lims
    
    lat_values = np.round(data.latitudes()[0][0::int((lon_max - lon_min)/reso)],1)
    lon_values = np.round(data.longitudes()[0][0:int((lon_max - lon_min)/reso)],1)

    lats_filtered_indices = [i for i in range(len(lat_values)) if lat_values[i] >= lat_low and lat_values[i] <= lat_high]
    lons_filtered_indices = [i for i in range(len(lon_values)) if lon_values[i] >= lon_low and lon_values[i] <= lon_high]
    
    return_data_struct = {}
    
    for z in levels[::-1]:
        uz = data.select(shortName = 'u',level=[z])
        vz = data.select(shortName = 'v',level=[z])
        
        u_filtered = []
        v_filtered = []

        grid =(int((lat_max - lat_min)/reso + 1), int((lon_max - lon_min)/reso))
        uz = reshape(uz.values(),grid)   
        vz = reshape(vz.values(),grid)
        
        u_lat_filter = uz[lats_filtered_indices[::-1]]
        v_lat_filter = vz[lats_filtered_indices[::-1]]
        
        for i in range(len(u_lat_filter)):
            u_filtered.append(u_lat_filter[i][lons_filtered_indices])
            v_filtered.append(v_lat_filter[i][lons_filtered_indices])
        
        return_data_struct[z] = (array(u_filtered),array(v_filtered))
        
    return return_data_struct, lat_values[lats_filtered_indices[::-1]],lon_values[lons_filtered_indices]


#linear interpolation on 3D grid
def interpol_linear(data, lat, lon, z, levels, lats, lons):
    
    u_grid_values = array([list(data.values())[i][0] for i in range(len(levels))])
    v_grid_values = array([list(data.values())[i][1] for i in range(len(levels))])

    u = interpn((levels[::-1],lats,lons),u_grid_values,[z,lat,lon])[0]
    v = interpn((levels[::-1],lats,lons),v_grid_values,[z,lat,lon])[0]

    return (u,v)
    


