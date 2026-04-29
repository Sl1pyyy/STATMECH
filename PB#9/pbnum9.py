import numpy as np
import matplotlib.pyplot as plt

#Problem 1
T = 3530 #K
R = 8.314 #J/(mol * K)

def fract(w, dE):
    """Nominator of fractions computation"""
    return w * np.exp(-dE /(R * T))

dE1 = 2740 #J/mol
dE2 = 19900 #J/mol
w0 = 9
w1 = 5
w2 = 7

G = fract(w0, 0)
F1 = fract(w1, dE1)
F2 = fract(w2, dE2)

#Partition function
Q = G + F1 + F2

#Fractions calculation
f1 = (F1 / Q) * 100
f2 = (F2 / Q) * 100
g = (G / Q) * 100
print(f'Fractions of Ce within temperature {T} K in first degenerate: {f1:.2f}%, second degenerate: {f2:.2f}% states'
      f'\nAt ground state fraction is = {g:.2f}%')

#Problem 2

#Electronic partition function
q = w0 + F1

#Total energy calculations

Ut = 1.5 * R * T
Uel = (dE1 * F1) / q

#Heat capacity contributions
Cvt = 1.5 * R

Cvel = R * (dE1 / (R * T)) * w0 * F1 / (q) ** 2

print(f'Electronic partition function = {q}, translational energy = {Ut:.2f} J/mol, electronic energy = {Uel:.2f} J/mol'
      f'\nHeat capacity (Cv,t) = {Cvt:.2f} J/(mol * K), heat capacity (Cv,el) = {Cvel:.2f} J/(mol * K)'
      f'\nRelative electronic contributions: energy - {Uel/(Uel  + Ut) * 100:.2f}%, heat capacity - {Cvel/(Cvel + Cvt) * 100:.2f}%')

#Problem 3
hbar = 1.0546 * 10 ** (-34) #J * s
kd = 320 #N/m
d0  = 1.99 * 10 ** (-10) #  m
T1 = 299 #K
T2 = 1000 #K
kB = 1.3806 * 10 ** (-23) # J/K
#Reduced mass Cl35Cl37
u = (35 * 37 / 72) * (1.66054 * 10 ** (-27)) #kg

#Vibrational temperature
wvib = (kd / u) ** 0.5

Ov = (hbar * wvib) / kB

#Rotational temperature
I = u * (d0 ** 2)

Or = (hbar ** 2) / (2 * I * kB)

print(f'\nVibrational temperature = {Ov:.2f} K, rotational temperature = {Or:.2f} K')

#Distribution plotting

def plot_distr(T):
    """Calculate vibrational distribution: P(v) = exp(-v * Ov / T) * (1 - exp(-Ov / T))
    and Rotational distribution: P(J) = (Or / T) * (2J + 1) * exp(-J(J+1) * Or / T)"""
    v = np.arange(0, 15)
    p_v = np.exp(-v * Ov / T) * (1 - np.exp(-Ov / T))
    cum_v = np.cumsum(p_v)

    j_limit = int(np.sqrt(T / Or) * 3)
    j = np.arange(0, j_limit)
    p_j = (Or / T) * (2*j + 1) * np.exp(-j * (j + 1) * (Or / T))
    cum_j = np.cumsum(p_j)


    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f'Distributions for 35Cl37Cl at T = {T} K', fontsize=16)

    # Vibrational plots
    axs[0, 0].bar(v, p_v, color='skyblue', edgecolor='black')
    axs[0, 0].set_title('Vibrational distribution P(v)')
    axs[0, 1].step(v, cum_v, where='mid', color='blue')
    axs[0, 1].set_title('Vibrational cumulative distribution')

    # Rotational plots
    axs[1, 0].plot(j, p_j, color='salmon')
    axs[1, 0].set_title('Rotational distribution P(J)')
    axs[1, 1].plot(j, cum_j, color='red')
    axs[1, 1].set_title('Rotational cumulative distribution')

    for ax in axs.flat:
        ax.set_xlabel('State number (v or J)')
        ax.grid(alpha=0.3)

    plt.show()

plot_distr(T1)
plot_distr(T2)

#Problem 4
D0 = 57.6 # kcal/mol
T4 = 800 #K

u2 = (35 / 2) * (1.66054 * 10 ** (-27)) #kg
wvib2 = (kd / u2) ** 0.5
Ov2 = hbar * wvib2 / kB
print(f'\nVibrational temperature = {Ov2:.2f} K of Cl2')

#Utotal = Utr + Urot + Uvib + Uel
Utotal = (1.5 * R * T4) + R * T4 + ((R * Ov2) / (np.exp(Ov / T4) - 1)) + 0
print(f'\nTotal energy contribution is = {Utotal:.2f} J/mol'
      f'\nAdding zero-point energy, total energy contribution is now = {Utotal + (R * Ov2) * 0.5:.2f} J/mol')

#Cv,vib = Ctr + Crot + Cvib + Cel
Cv2 = 1.5 * R + R + R * ((Ov2 / T4) ** 2) * (np.exp(Ov2 / T4) / (np.exp(Ov2 / T4) - 1) ** 2)
print(f'\nTotal heat capacity contribution is = {Cv2:.2f} J/(mol * K)')






