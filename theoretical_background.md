# Trajectory simulation
## Ascent model
$$m_{tot}a_z = g\rho_{air}V - g(m_{gross} + m_{gas}) - \frac{1}{2}C_D\rho_{air}\dot{z}^2A_b$$
where:\
$m_{tot} = m_{gross} + m_{gas} + m_{added}$\
$m_{added} = \frac{1}{2}\rho_{air}V$ - the mass of the air
moved by the balloon during the ascent\
**The rise in volume is the main responsible for the rise in
the balloon vertical velocity**, commonly called ascending
rate, because in the equation of motion the acceleration is observed to depend on volume.

$$V = \frac{m_{gas}}{M_{gas}} R \frac{T_{gas}}{p_{air}}$$

Since our approximation addresses the volume change totally to the adiabatic expansion:

$$p \cdot V^{\gamma} = constant$$


### Other related stuff
$C_D$ - so-called **drag coefficient** $\rightarrow$ yet to be figured out (or guesstimated ;) )\
**sources**: https://www.researchgate.net/publication/253549900_Modeling_the_ascent_of_sounding_balloons_Derivation_of_the_vertical_air_motion \
$A_b$ - balloon reference area ($m^2$)

#### Altitude variation of the pressure:
Theoretically valid only up to $\sim$ 20km

$$p = p_0 \cdot \left(1 - \frac{L \cdot h}{T_0} \right)^{\frac{g \cdot M}{R_0 \cdot L}} = p_0 \cdot \left(1 - \frac{g \cdot h}{c_p \cdot T_0} \right)^{\frac{c_p \cdot M}{R_0}} \approx p_0 \cdot exp \left(-\frac{g \cdot h \cdot M}{T_0 \cdot R_0} \right)$$
Check out the following link if you are interested in what the heck those parameters represent (actually they are straightforward...):\
https://en.wikipedia.org/wiki/Atmospheric_pressure?fbclid=IwAR2keT2FBSF0ExzFWVDx_fpQHDhbOnQdD49PGO7yijUqn8QY_0d6BH12uc0

#### Altitude variation of the density:
Only up to $\sim$ 20 km? It maybe won`t be enough. **We should discuss it next time!**

$$\rho = \frac{p \cdot M}{R \cdot T}$$
where:\
$p$ is described previously\
$T = T_0 - L \cdot h$

$$\rho = \frac{p_0 \cdot M}{R \cdot T_0} \left( 1 - \frac{L \cdot h}{T_0} \right)^{\frac{g \cdot M}{R \cdot L} - 1}$$

**Source:** https://en.wikipedia.org/wiki/Density_of_air

## Descent model
After the balloon bursts, the payload starts to drop slowed
down by a parachute. The system acts as it is not subjected
to inertial acceleration: after a brief transient of time weight
is perfectly balanced by aerodynamic drag; then, the
payload will fall down at an approximately constant speed
called terminal velocity

$$v_t = \sqrt{\frac{2 \cdot m \cdot g}{C_D \cdot \rho_{air} \cdot A_p}}$$

## Laterla displacement
"As the wind blows..." and also probably the hardest to figure out...
