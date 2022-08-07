import numpy as np
import matplotlib.pyplot as plt
from wind_velocity import *
import pandas
import sys 

#wind data for horizontal displacement
data = ecdata.read("20220805000000-6h-oper-fc.grib2")
levels = [1000, 925, 850, 700, 500, 300, 250, 200, 50] #altitude in mbars
lat_lims = (40,50) #latitude limits
lon_lims = (20,30) #longitud limits

dataf,lats,lons = filter_to_region(data,lat_lims,lon_lims,levels) #filtering according to the previously given limits

#46.7500, 22.9667 Meregyo apporximate coordinates of the launching site

z_0 = 0. #initial altitude
v_0 = 0. #initial velocity
lat_0 = 46.7500 #initial latitude
lon_0 = 22.9667 #initial longitude

R_earth = 6317000 #earth mean radius
g = 9.81 #gravitational acceleration
C_d_balloon = 0.25 #drag coefficient of the balloon
C_d_parachute = 1.75 #drag coefficient of the parachute
A_p = np.pi/4 #parachute reference area
V_max = 33. #maximal balloon volume (m^3)
V_0 = 2.3 #initial volume of the balloon (m^3)
rho_helium = 0.1785 #(kg/m^3)

#masses
m_gross = 1.5 #payload mass + balloon mass
m_gas = V_0*rho_helium #mass of the helium
m_added = 0. #mass of the additional air, it`s value is 0, since we use pure helium

m_tot = m_gross + m_gas + m_added #total mass

dt = 1. #time step in seconds

M_air = 28.97*10**(-3) #mean molar mass of dry (kg/mol) I have to express it this way to obtain the same dimensions everywhere
R = 8.314 #universal gas constant (J/(K*mol))
R_air = R/M_air #gas constant for air

T_0 = 288.16 #sea level standard temperature (K)
p_0 = 101325 #sea level standard atmospheric pressure (Pa) N/m^2
rho_0 = 1.225 #density of air at sea level (kg/m^3)

#boundary heights in m according to the "US Standard Atmosphere 1976": https://web.archive.org/web/20070310223946/http://www.atmosculator.com/The%20Standard%20Atmosphere.html
h_0 = 0.
h_1 = 11000.
h_2 = 20000.
h_3 = 32000.

#lapse rates in the different domains in (K/m)
lambd_0 = -6.5/1000 
lambd_1 = 0./1000
lambd_2 = 1.0/1000

def temperature(z_n):
	if((z_n >= 0) and (z_n <= h_1)):
		return T_0 + z_n*lambd_0
	if((z_n > h_1) and (z_n <= h_2)):
		return temperature(h_1)
	if((z_n > h_2) and (z_n <= h_3)):
		return temperature(h_2) + (z_n - h_2)*lambd_2

def pressure(z_n):
	if((z_n >= 0) and (z_n <= h_1)):
		return p_0*(1 + z_n*lambd_0/T_0)**(-g/lambd_0/R_air)
	if((z_n > h_1) and (z_n <= h_2)):
		return pressure(h_1)*np.exp(-(z_n - h_1)*g/(R_air*temperature(h_1)))
	if((z_n > h_2) and (z_n <= h_3)):
		return pressure(h_2)*(temperature(z_n)/temperature(h_2))**(-g/(lambd_2*R_air))

def rho_air(z_n):
	if((z_n >= 0) and (z_n <= h_1)):
		return rho_0*(1 + z_n*lambd_0/T_0)**(-g/lambd_0/R_air - 1)
	if((z_n > h_1) and (z_n <= h_2)):
		return rho_0*(pressure(h_1)*T_0)/(p_0*temperature(h_1))*np.exp(-(z_n - h_1)*g/R_air/temperature(h_1))
	if((z_n > h_2) and (z_n <= h_3)):
		return (rho_0/p_0)*pressure(h_2)*( temperature(h_2)/T_0 + (z_n - h_2)*lambd_2/T_0 )**(-g/(lambd_2*R_air) - 1) * (T_0/temperature(h_2))**(-g/(lambd_2*R_air))


# volume from the ideal gas law assuming PERFECT thermal equilibrium with the atmosphere
# V = Nu*R*T/p

mu_helium = 4 #molar mass of Helium (kg/kmol)
Nu = m_gas/mu_helium*1000 #mol
#print("V_0_ideal_gas = ", Nu*R*T_0/p_0)

def volume(z_n):
	return Nu*R*temperature(z_n)/pressure(z_n)

'''
# volume assuming an adiabatic process
gamma = 5/3 #adiabatic coefficient for helium (ideal monoatomic gas!!!)
const = pressure(z_0)*V_0**gamma #p*V^gamma = const

def volume(z_n):
	return (const/pressure(z_n))**(1/gamma)
'''
'''
# volume assuming ploytropic process according to: https://upcommons.upc.edu/bitstream/handle/2117/183322/Report_fitxer%20de%20consulta.pdf?fbclid=IwAR3G-xUTHAHRYchAfr6tSpMtpGVDH9EVD0C8rjRiJsCGgepaciWBHw-MNl8
def volume(z_n):
	if(z_n <= altitude_c):
		return (const1/pressure(z_n))**(1/n1)
	else:
		return (const2/pressure(z_n))**(1/n2)

altitude_c = 9800 #poly1 to poly2 domain edge (or whatever... you know what I mean)
n1 = 0.7 #polytropic coefficient in the first domain
const1 = pressure(z_0)*V_0**n1
n2 = 0.9 #polytropic coefficient in the second domain
const2 = pressure(altitude_c)*volume(altitude_c)**n2
'''
def A_b_balloon(z_n):
	return np.pi*((3/(4*np.pi)*volume(z_n))**(1/3))**2

################################################################ Ascent ########################################################################

f_1 = lambda z_n, v_n : (1/m_tot)*(g*rho_air(z_n)*volume(z_n) - g*(m_gross + m_gas) - 1/2*(C_d_balloon*rho_air(z_n)*v_n**2*A_b_balloon(z_n)))

#there is no explicit time dependence
def velocity_ascending(z_n, v_n):
	K_1 = f_1(z_n, v_n)
	K_2 = f_1(z_n, v_n + dt*K_1/2)
	K_3 = f_1(z_n, v_n + dt*K_2/2)
	K_4 = f_1(z_n, v_n + dt*K_3)
	return v_n + 1/6*(K_1 + 2*K_2 + 2*K_3 + K_4)*dt

################################################################ Descent #######################################################################

def velocity_descending(z_n, v_n):
	return -np.sqrt(2*m_gross*g/(C_d_parachute*rho_air(z_n)*A_p))

### Altitude ###

dz = lambda dt, v_n : dt*v_n 

f_2 = lambda z_n, v_n : v_n

#there is no explicit time dependence
def altitude(z_n, v_n):
	K_1 = f_2(z_n, v_n)
	K_2 = f_2(z_n + dz(dt, v_n)*K_1/2, v_n)
	K_3 = f_2(z_n + dz(dt, v_n)*K_2/2, v_n)
	K_4 = f_2(z_n + dz(dt, v_n)*K_3, v_n)
	return z_n + 1/6*(K_1 + 2*K_2 + 2*K_3 + K_4)*dt

#initial conditions
z_n = z_0
v_n = v_0
lat_n = lat_0
lon_n = lon_0

V = np.array([v_0]) #array for the velocities
Z = np.array([z_0]) #array for the altitudes
Lat = np.array([lat_0]) #array for the latitudes
Lon = np.array([lon_0]) #array for the longitudes
U_wind = np.array([]) #array for the "u" wind components
V_wind = np.array([]) #array for the "v" wind components
Volume = volume(0.)
Coordinates = np.array([]) #array for the coordinates, example for a row:(lon, lat, z)

######## Ascending phase ########
while(Volume < V_max): #iterates up to the burst volume of the balloon
	print("z_n = ", z_n)
	v_np1 = velocity_ascending(z_n, v_n) #calculating the new velocity
	z_np1 = altitude(z_n, v_n) #calculating the new altitude
	v_n = v_np1
	z_n = z_np1
	z_n_mbar = pressure(z_n)/100 #in order to get mbar, to use the wind data

	u_wind, v_wind = interpol_linear(dataf, lat_n, lon_n, z_n_mbar, levels, lats, lons) #extracting the wind components
	U_wind = np.append(U_wind, u_wind)
	V_wind = np.append(V_wind, v_wind)
	
	lat_n += (180/np.pi)*v_wind*dt/R_earth #latitude displacement
	lon_n += (180/np.pi)*u_wind*dt/(R_earth*np.cos(lat_n*np.pi/180)) #longitude displacement

	Volume = volume(z_n)
	V = np.append(V, v_n)
	Z = np.append(Z, z_n)
	Lat = np.append(Lat, lat_n)
	Lon = np.append(Lon, lon_n)

######### Descending phase ########
while(z_n > z_0 ): #iterates while it reaches the surface
	print("z_n = ", z_n)
	v_np1 = velocity_descending(z_n, v_n)
	z_np1 = altitude(z_n, v_n)
	v_n = v_np1
	z_n = z_np1

	#if the integration goes under 0 altitude it is considered to be 0
	if(z_n < 0):
		z_n = 0
	
	z_n_mbar = pressure(z_n)/100 #in order to get mbar, to use the wind data

	u_wind, v_wind = interpol_linear(dataf, lat_n, lon_n, z_n_mbar, levels, lats, lons) #extracting the wind components
	U_wind = np.append(U_wind, u_wind)
	V_wind = np.append(V_wind, v_wind)

	lat_n += (180/np.pi)*v_wind*dt/R_earth #latitude displacement
	lon_n += (180/np.pi)*u_wind*dt/(R_earth*np.cos(lat_n*np.pi/180)) #longitude displacement

	Volume = volume(z_n)
	V = np.append(V, v_n)
	Z = np.append(Z, z_n)
	Lat = np.append(Lat, lat_n)
	Lon = np.append(Lon, lon_n)

#creating the coordinate data structure
Coordinates = np.concatenate((Lat.reshape(-1,1), Lon.reshape(-1,1), Z.reshape(-1,1)), axis = 1) 
print(Coordinates)
pandas.DataFrame(data = Coordinates, columns = ["latitude", "longitude", "altitude"]).to_csv("Coordinates_test.csv", index = False)

#### Representing the altitude time series
plt.plot(Z)
plt.grid()
plt.xlabel('time (s)')
plt.ylabel('altitude (m)')
plt.title('Altitude in time')
plt.show()