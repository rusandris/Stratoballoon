import numpy as np
import matplotlib.pyplot as plt

sim_length = 1000

z_0 = 0. #initial altitude
v_0 = 0. #initial velocity

#masses
m_gross = 1. 
m_gas = 1.
m_added = 1.

m_tot = m_gross + m_gas + m_added

rho_air = 1. #air density

g = 10. #gravitational acceleration
C_d = 1. #drag coefficient
A_b = 1. #balloon reference area

dt = 1. #time step

M_gas = 4*10**(-3) #molar mass (kg/mol) I have to express it this way to obtain the same dimensions everywhere
R = 8.314 #universal gas constant (J/(K*mol))

#pressure parameters
c_p = 1004.68506 #constant pressure specific heat (J/(kg/K))
T_0 = 288.16 #sea level standard temperature (kg)
M = 0.02896968 #molar mass of dry air (kg/mol)
p_0 = 101325 #sea level standard atmospheric pressure (Pa)  N/m^2

def pressure(z_n):
	return p_0*(1 - (g*z_n)/(c_p*T_0))**(c_p*M/R)


V_0 = 2.34 #initial volume of the balloon (m^3)
gamma = 5/3 #adiabatic coefficient for helium (ideal monoatomic gas!!!)

const = p_0*V_0**gamma #p*V^gamma = const

def volume(z_n):
	return (const/pressure(z_n))**(1/gamma)

def velocity(z_n, v_n):
	#print(v_n)
	return (1/m_tot)*(m_tot*v_n + g*rho_air*volume(z_n)*dt - g*(m_gross + m_gas)*dt - 1/2*(C_d*rho_air*v_n**2*A_b*dt))

def altitude(z_n, v_n):
	return z_n + v_n*dt

z_n = z_0
v_n = v_0

V = np.array([])
Z = np.array([])
for i in range(sim_length):
	V = np.append(V, v_n)
	Z = np.append(Z, z_n)
	v_np1 = velocity(z_n, v_n)
	z_np1 = altitude(z_n, v_n)
	v_n = v_np1
	z_n = z_np1

print(V)
print(Z)

plt.plot(Z)
plt.xlabel('time...sort of...')
plt.ylabel('altitude...kindof...')
plt.show()