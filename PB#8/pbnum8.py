import matplotlib.pyplot as plt
import numpy as np
from openpyxl.chart import label

h = 6.626 * 10 ** (-34)  #J*s
kB = 1.3806 * 10 ** (-23)  # J/K

#Source of experimental data: C. Kittel, Introduction to Solid State Physics (2005)
#For the velocities were used Debye average velocity, which is the sum of longitudinal and traverse velocities,
#due to  the fact, that using only tabular longitudinal velocity will lead to overestimation of Debye temperature
v0c = 2611 # m/s
v0lc = 4760 #m/s
v0d = 13950 # m/s
v0ld = 17500 #m/s
roc = (8960 * 6.022 * 10 ** 23) / 0.063546 # atoms/m^3
rod = (3512 * 6.022 * 10 ** 23) / 0.012011 # atoms/m^3

#Problem 1
def debtemp (v0, ro):
    return (h * v0 * (6 * (np.pi) ** 2 * ro) ** (1/3) )/ (2 * np.pi * kB)

copper = debtemp(v0c, roc)
copper_long = debtemp(v0lc, roc)
print(f'Copper Debye temperature is equal to {copper:.2f} K'
      f'\nSource says that the copper Debye temperature is equal to 343 K, which is very close to calculated one. '
      f'Using only longitudinal velocity will give: {copper_long:.2f} K ')

diamond = debtemp(v0d, rod)
diamond_long = debtemp(v0ld, rod)
print(f'Diamond Debye temperature is equal to {diamond:.2f} K'
      f'\nSource says that the diamond Debye temperature range is 2240 K, hence calculated  one is overestimated by {2240 - diamond:.2f} K. '
      f'Using only longitudinal velocity will give:  {diamond_long:.2f} K ')


#Problem 3
r_values = np.linspace(0.001, 2.0, 100)
t = np.linspace(0.001, 2.0, 100)

from scipy.integrate import quad

#Define the integrand function
def integrand(x):
    return (x**4 * np.exp(x)) / (np.exp(x) - 1)**2

#Define the Cv function and calculate it using Gaussian quadrature
def debye_cv(r_ratio):
    upper_limit = 1.0 / r_ratio
    integral, _ = quad(integrand, 0, upper_limit)
    return 9 * (r_ratio**3) * integral

cv_values = [debye_cv(r) for r in r_values]
debye = t ** 3  * ((12 * np.pi ** 4) / 5 )

plt.plot(r_values, cv_values, label = ' Heat capacity in T/ΘD')
plt.plot(r_values, debye, label = 'Debye T^3 law')
plt.axhline(y = 3, color = 'purple', linestyle = '--', label = 'Dulong-Petit limit 3')
plt.ylim(top = 3.1, bottom = 0)
plt.xlabel('T/ΘD')
plt.ylabel('Cv/R')
plt.legend()
plt.grid(True)
plt.show()



