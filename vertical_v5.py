import numpy as np
import matplotlib.pyplot as plt
from wind_velocity import *
import pandas
import sys 

print(sys.argv)

#wind data for horizontal displacement
data = ecdata.read("20220805000000-6h-oper-fc.grib2")
levels = [1000, 925, 850, 700, 500, 300, 250, 200, 50]
lat_lims = (40,50)
lon_lims = (20,30)

dataf,lats,lons = filter_to_region(data,lat_lims,lon_lims,levels)

#46.590743, 22.923171
#46.7500, 22.9667 Meregyo

z_0 = 0. #initial altitude
v_0 = 0. #initial velocity
lat_0 = 46.7500 #initial latitude
lon_0 = 22.9667 #initial longitude

R_earth = 6317000 #earth mean radius

g = 9.81 #gravitational acceleration
C_d_balloon = 0.25 #drag coefficient of the balloon
C_d_parachute = 1.75 #drag coefficient of the parachute
A_b = 1. #balloon reference area, changes in time!!!!
A_p = np.pi/4 #parachute reference area
V_max = 30. #maximal balloon volume (m^3)
V_0 = 2.3 #initial volume of the balloon (m^3)
rho_helium = 0.169 #(kg/m^3)
#masses
m_gross = 1.5
m_gas = V_0*rho_helium
m_added = 0.

m_tot = m_gross + m_gas + m_added

dt = 1.

M_air = 28.97*10**(-3) #mean molar mass of dry (kg/mol) I have to express it this way to obtain the same dimensions everywhere
R = 8.314 #universal gas constant (J/(K*mol))
R_air = R/M_air

T_0 = 288.16 #sea level standard temperature (K)
p_0 = 101325 #sea level standard atmospheric pressure (Pa)  N/m^2
rho_0 = 1.225 #density of air at sea level (kg/m^3)

#boundary heights in m
h_0 = 0.
h_1 = 11000.
h_2 = 20000.
h_3 = 32000.

#lapse rates in the different domains
lambd_0 = -6.5/1000 #(J/m)!!!!
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
	#print(z_n)
	if((z_n >= 0) and (z_n <= h_1)):
		return rho_0*(1 + z_n*lambd_0/T_0)**(-g/lambd_0/R_air - 1)
	if((z_n > h_1) and (z_n <= h_2)):
		return rho_0*(pressure(h_1)*T_0)/(p_0*temperature(h_1))*np.exp(-(z_n - h_1)*g/R_air/temperature(h_1))
	if((z_n > h_2) and (z_n <= h_3)):
		return (rho_0/p_0)*pressure(h_2)*( temperature(h_2)/T_0 + (z_n - h_2)*lambd_2/T_0 )**(-g/(lambd_2*R_air) - 1) * (T_0/temperature(h_2))**(-g/(lambd_2*R_air))


# volume from the ideal gas law
# V = Nu*R*T/p

mu_helium = 4 #molar mass of Helium (kg/kmol)
Nu = m_gas/mu_helium*1000 #mol
#print("V_0_ideal_gas = ", Nu*R*T_0/p_0)

def volume(z_n):
	return Nu*R*temperature(z_n)/pressure(z_n)

'''
# volume assuming an adiabatic process
gamma = 5/3 #adiabatic coefficient for helium (ideal monoatomic gas!!!)

const = p_0*V_0**gamma #p*V^gamma = const

def volume(z_n):
	return (const/pressure(z_n))**(1/gamma)
'''
def A_b_balloon(z_n):
	return (3/(4*np.pi)*volume(z_n))**(1/3)
################################################################ Ascent ########################################################################
'''
def velocity_ascending(z_n, v_n):
	print('v_n = ', v_n)
	print('pressure = ', pressure(z_n))
	print('volume = ', volume(z_n))
	return (1/m_tot)*(m_tot*v_n + g*rho_air(z_n)*volume(z_n)*dt - g*(m_gross + m_gas)*dt - 1/2*(C_d*rho_air(z_n)*v_n**2*A_b*dt))
'''
#implemeting the same with RK4 method
dz = lambda dt, v_n : dt*v_n 

f = lambda z_n, v_n : (1/m_tot)*(g*rho_air(z_n)*volume(z_n) - g*(m_gross + m_gas) - 1/2*(C_d_balloon*rho_air(z_n)*v_n**2*A_b_balloon(z_n)))

#there is no explicit time dependence
def velocity_ascending(z_n, v_n):
	K_1 = f(z_n, v_n)
	#print("K_1 = ", K_1)
	K_2 = f(z_n, v_n + dt*K_1/2)
	K_3 = f(z_n, v_n + dt*K_2/2)
	K_4 = f(z_n, v_n + dt*K_3)
	return v_n + 1/6*(K_1 + 2*K_2 + 2*K_3 + K_4)*dt

################################################################ Descent #######################################################################

def velocity_descending(z_n, v_n):
	return -np.sqrt(2*m_gross*g/(C_d_parachute*rho_air(z_n)*A_p))

'''
def altitude(z_n, v_n):
	return z_n + v_n*dt
'''

f_2 = lambda z_n, v_n : v_n

#implementing the same with RK4 method
#ez itt még hibás!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def altitude(z_n, v_n):
	K_1 = f_2(z_n, v_n)
	#print("K_1 = ", K_1)
	K_2 = f_2(z_n + dz(dt, v_n)/2, v_n)
	K_3 = f_2(z_n + dz(dt, v_n)/2, v_n)
	K_4 = f_2(z_n + dz(dt, v_n), v_n)
	return z_n + 1/6*(K_1 + 2*K_2 + 2*K_3 + K_4)*dt


z_n = z_0
v_n = v_0
lat_n = lat_0
lon_n = lon_0

V = np.array([v_0])
Z = np.array([z_0])
Lat = np.array([lat_0])
Lon = np.array([lon_0])
U_wind = np.array([])
V_wind = np.array([])
Volume = volume(0.)
Coordinates = np.array([])

while(Volume < V_max): 
	'''print('z_n kint = ', z_n)
	print('v_n kint = ', v_n)
	print('V kint = ', Volume)
	'''
	V = np.append(V, v_n)
	Z = np.append(Z, z_n)
	v_np1 = velocity_ascending(z_n, v_n)
	z_np1 = altitude(z_n, v_n)
	v_n = v_np1
	z_n = z_np1
	z_n_mbar = pressure(z_n)/100 #(in order to get mbar)
	#print('z_n_mbar = ', z_n_mbar)

	u_wind, v_wind = interpol_linear(dataf, lat_n, lon_n, z_n_mbar, levels, lats, lons)
	U_wind = np.append(U_wind, u_wind)
	V_wind = np.append(V_wind, v_wind)
	
	lat_n += (180/np.pi)*v_wind*dt/R_earth
	lon_n += (180/np.pi)*u_wind*dt/(R_earth*np.cos(lat_n*np.pi/180))

	Volume = volume(z_n)
	Lat = np.append(Lat, lat_n)
	Lon = np.append(Lon, lon_n)
	'''
	print("lat = ", lat_n)
	print("lon = ", lon_n)
	'''
while(z_n > z_0 + 10): #clarify it!!!!!!!! #+10 not to reach negative values

	v_np1 = velocity_descending(z_n, v_n)
	z_np1 = altitude(z_n, v_n)
	v_n = v_np1
	z_n = z_np1
	'''
	print('z_n kint = ', z_n)
	print('v_n kint = ', v_n)
	print('V kint = ', Volume)
	'''
	z_n_mbar = pressure(z_n)/100 #(in order to get mbar)
	'''
	print('z_n_mbar_test = ', pressure(z_n)/100)
	print('z_n_mbar = ', z_n_mbar)
	'''
	u_wind, v_wind = interpol_linear(dataf, lat_n, lon_n, z_n_mbar, levels, lats, lons)
	U_wind = np.append(U_wind, u_wind)
	V_wind = np.append(V_wind, v_wind)

	lat_n += (180/np.pi)*v_wind*dt/R_earth
	lon_n += (180/np.pi)*u_wind*dt/(R_earth*np.cos(lat_n*np.pi/180))
	Volume = volume(z_n)
	V = np.append(V, v_n)
	Z = np.append(Z, z_n)
	Lat = np.append(Lat, lat_n)
	Lon = np.append(Lon, lon_n)
	'''
	print("lat = ", lat_n)
	print("lon = ", lon_n)
	'''
print(V)
print(Z)

Coordinates = np.concatenate((Lat.reshape(-1,1), Lon.reshape(-1,1), Z.reshape(-1,1)), axis = 1)
print(Coordinates)

pandas.DataFrame(data = Coordinates, columns = ["latitude", "longitude", "altitude"]).to_csv("Coordinates_test.csv", index = False)

plt.plot(Z)
plt.grid()
plt.xlabel('time...sort of...')
plt.ylabel('altitude')
plt.title('Altitude in time')
plt.show()