import numpy as np
import matplotlib.pyplot as plt

z_0 = 0. #initial altitude
v_0 = 0. #initial velocity

#masses
m_gross = 1. 
m_gas = 1.
m_added = 1.

m_tot = m_gross + m_gas + m_added

g = 9.81 #gravitational acceleration
C_d = 1. #drag coefficient
A_b = 1. #balloon reference area, changes in time!!!!
A_p = np.pi/4 #parachute reference area
V_max = 30 #maximal balloon volume (m^3)

dt = 0.01

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

'''
# volume from the ideal gas law
# V = Nu*R*T/p

Nu = 89.2 #(mol)

def volume(z_n):
	return Nu*R*temperature(z_n)/pressure(z_n)
'''

# volume assuming an adiabatic process
V_0 = 2.34 #initial volume of the balloon (m^3)
gamma = 5/3 #adiabatic coefficient for helium (ideal monoatomic gas!!!)

const = p_0*V_0**gamma #p*V^gamma = const

def volume(z_n):
	return (const/pressure(z_n))**(1/gamma)


################################################################ Ascent ########################################################################

def velocity_ascending(z_n, v_n):
	print('v_n = ', v_n)
	print('pressure = ', pressure(z_n))
	print('volume = ', volume(z_n))
	return (1/m_tot)*(m_tot*v_n + g*rho_air(z_n)*volume(z_n)*dt - g*(m_gross + m_gas)*dt - 1/2*(C_d*rho_air(z_n)*v_n**2*A_b*dt))

################################################################ Descent #######################################################################

def velocity_descending(z_n, v_n):
	return -np.sqrt(2*m_gross*g/(C_d*rho_air(z_n)*A_p))

def altitude(z_n, v_n):
	return z_n + v_n*dt


z_n = z_0
v_n = v_0

V = np.array([])
Z = np.array([])
Volume = volume(0.)
while(Volume < V_max): 
	print('z_n kint = ', z_n)
	V = np.append(V, v_n)
	Z = np.append(Z, z_n)
	v_np1 = velocity_ascending(z_n, v_n)
	z_np1 = altitude(z_n, v_n)
	v_n = v_np1
	z_n = z_np1
	'''
	x_wind = interpol_x(x_n, y_n, z_n)
	x_n += x_wind*dt/R
	y_n += 
	'''
	Volume = volume(z_n)

while(z_n > 0):
	print('z_n kint = ', z_n)
	V = np.append(V, v_n)
	Z = np.append(Z, z_n)
	v_np1 = velocity_descending(z_n, v_n)
	z_np1 = altitude(z_n, v_n)
	v_n = v_np1
	z_n = z_np1

print(V)
print(Z)

plt.plot(Z)
plt.grid()
plt.xlabel('time...sort of...')
plt.ylabel('altitude')
plt.title('Altitude in time')
plt.show()