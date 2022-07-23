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

### Other related stuff
$C_D$ - so-called **drag coefficient** $\rightarrow$ yet to be figured out\
**sources**: https://www.researchgate.net/publication/253549900_Modeling_the_ascent_of_sounding_balloons_Derivation_of_the_vertical_air_motion \
$A_b$ - balloon reference area ($m^2$)

#### Altitude variation of the pressure:

$$p = p_0 \cdot(1 - \frac{g \cdot h}{c_p \cdot T_0})^{\frac{c_p \cdot M}{R_0}} \approx p_0 \cdot exp \left(-\frac{g \cdot h \cdot M}{T_0 \cdot R_0} \right)$$
