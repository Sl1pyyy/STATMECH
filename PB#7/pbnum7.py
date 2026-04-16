import numpy as np
import matplotlib.pyplot as plt

h = 6.626e-34  #J*s
c = 2.9979e8  #m/s
kB = 1.3806e-23  # J/K


def planck(v, T):
    """Calculates spectral energy density u_v (J*s/m^3)"""
    term1 = (8 * np.pi * h * v ** 3) / (c ** 3)
    exponent = (h * v) / (kB * T)
    return term1 * (1 / (np.exp(exponent) - 1))

temps = [373, 1373, 1000000]
labels = ['Boiling Water (373K)', 'Burning Charcoal (1373K)', 'Sun Corona (10^6K)']
colors = ['royalblue', 'darkorange', 'crimson']

#Define Frequency Range (Infrared to X-ray)
v = np.logspace(11, 19, 1000)

plt.figure(figsize=(12, 8))

for T, label, col in zip(temps, labels, colors):
    u_v = planck(v, T)
    plt.plot(v, u_v, label=label, color=col, lw=2.5)

    #Calculate and Point Out Most Probable Frequency (Wien's Law for Frequency)
    v_max = (2.821 * kB * T) / h
    u_max = planck(v_max, T)

    # Mark the peak with a dot
    plt.scatter(v_max, u_max, color='black', zorder=5)

    # Annotate the peak frequency
    plt.annotate(f'Peak: {v_max:.2e} Hz',
                 xy=(v_max, u_max),
                 xytext=(0, 15),
                 textcoords='offset points',
                 ha='center',
                 fontsize=10,
                 fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.7))

plt.xscale('log')
plt.yscale('log')
plt.title("Planck's Law: Spectral Energy Density vs. Frequency", fontsize=14)
plt.xlabel('Frequency (Hz)', fontsize=12)
plt.ylabel('Energy density (J·s/m³)', fontsize=12)
plt.ylim(1e-25, 1e10)
plt.xlim(1e11, 1e19)
plt.legend(loc='upper left', fontsize=11)
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.tight_layout()
plt.show()