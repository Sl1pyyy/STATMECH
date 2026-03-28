import math
import numpy as np
import matplotlib.pyplot as plt

#Problem 1
g0, gp = 2, 1
T = 80 + 273.15
dE = 1.06
R = 1.987 * 10**-3

#n0/nt = g0 * exp(-E0/RT)/ exp(-E0/RT + -Ep/RT)
#np/nt = gp * exp(-Ep/RT)/ exp(-E0/RT + -Ep/RT)
#Dividied fist by second and got the resulting equation:
# n0/np = g0/gp * exp(-dE/RT)

ratio = (g0 / gp) * math.exp(-dE/(T*R))
print(f'The ratio is {ratio:.4f}')

#Fractions sum is equal 1: x0 + xp = 1
xp = 1/ (1 + ratio)
x0 = 1 - xp

print(f'The ortho is {x0:.4f} and para is {xp:.4f}')

#Problem 2
def distr(v, T, m):
    return (m * v / (kb * T)) * np.exp(-(m * v**2) / (2 * kb * T))

m = 1.66 * 10**(-27)
Rg =  8.3145
v = np.linspace(0, 7500, 500)
kb =  1.3803 * 10**(-23)

plt.plot(v, distr(v, 200, m), color='red', label='Distribution plot for 200K')
plt.plot(v, distr(v, 300, m), color='black', label='Distribution plot for 300K')
plt.xlabel('Velocity [m/s]')
plt.ylabel('Probability density function [s/m]')
plt.title('Boltzmann 2D distribution')
plt.grid(True)
plt.legend()
plt.show()
