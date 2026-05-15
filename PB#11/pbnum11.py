import numpy as np
import matplotlib.pyplot as plt
import math

R = 8.314 * 10 ** (-3)  # kJ / (mol * K)
H_cm = 1.4388 # K

def vib_energy(wavenumber, T):
    """Calculating vibrational internal energy contribution in kJ/mol"""
    e = 0
    for w in wavenumber:
        Ov = w * H_cm
        e += R * (Ov / (np.exp(Ov / T) -1))
        return e

#Problem 1

def xef2_reaction():
    """Xe + F2 = XeF2 reaction energy plot"""
    de_f2 = 154.8
    de_xef2 = 267.8
    w_f2 = [894]
    w_xef2 = [214, 214, 515, 555]

    T = np.linspace(300, 2000, 100)
    delta_u = []

    delta_elec = -(de_xef2 - de_f2)

    for t in T:
        u_trans = -1.5 * R * t
        u_vib = vib_energy(w_xef2, t) - vib_energy(w_f2, t)
        #delta_rot is 0 as derived in report
        delta_u.append(delta_elec + u_trans + u_vib)

    plt.figure(figsize=(10, 8))
    plt.plot(T, delta_u,color='skyblue')
    plt.title("Reaction energy for Xe + F2 reaction")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Energy change (U) [kJ/mol]")
    plt.grid(True)
    plt.show()

#Problem 2
def nitrogen_exchange():
    """N2 isotope exchange and Arrhenius plot."""
    O_14_14 = 3374.0

    m14 = 14.003
    m15 = 15.000

    mu_14_14 = (m14 * m14) / (m14 + m14)
    mu_15_15 = (m15 * m15) / (m15 + m15)
    mu_14_15 = (m14 * m15) / (m14 + m15)

    O_15_15 = O_14_14 * math.sqrt(mu_14_14 / mu_15_15)
    O_14_15 = O_14_14 * math.sqrt(mu_14_14 / mu_14_15)

    inv_temp_list = []
    log_k_list = []


    for T in range(300, 2050, 50):
        #Symmetry factor is (2*2)/(1*1) = 4
        symmetry = 4.0

        M1414 = 2 * m14
        M1515 = 2 * m15
        M1415 = m14 + m15
        mass_factor = ((M1415 ** 2) / (M1414 * M1515)) ** 1.5

        q1414 = math.exp(-O_14_14 / (2 * T)) / (1 - math.exp(-O_14_14 / T))
        q1515 = math.exp(-O_15_15 / (2 * T)) / (1 - math.exp(-O_15_15 / T))
        q1415 = math.exp(-O_14_15 / (2 * T)) / (1 - math.exp(-O_14_15 / T))

        vib_factor = (q1415 ** 2) / (q1414 * q1515)

        k_val = symmetry * mass_factor * vib_factor

        inv_temp_list.append(1.0 / T)
        log_k_list.append(math.log10(k_val))

    plt.plot(inv_temp_list, log_k_list, color="blue")
    plt.title("Arrhenius plot for nitrogen isotope exchange")
    plt.xlabel("1/T [K^-1]")
    plt.ylabel("log10(K)")
    plt.grid(True)
    plt.show()


#Problem 3

def N14_N15_exchange():
    # Each column is [Temperature, Time, N14N15_units, N15N15_units]
    data = [
        [465, 2, 100.5, 51.50],
        [465, 17, 125.8, 39.00],
        [465, 37, 163.0, 22.60],
        [465, 60, 194.0, 14.00],
        [500, 2, 108.0, 48.00],
        [500, 46, 207.0, 13.25]
    ]

    print("Temp(C) | Time(h) | Frac_14-14 | Frac_14-15 | Frac_15-15 | K_app")
    print("------------------------------------------------------------------")

    for row in data:
        temp = row[0]
        time = row[1]
        n1415 = row[2]
        n1515 = row[3]
        n1414 = 1000.0  #Reference base

        total_units = n1414 + n1415 + n1515
        f1414 = n1414 / total_units
        f1415 = n1415 / total_units
        f1515 = n1515 / total_units

        # Formula: K = [14-15]^2 / ([14-14] * [15-15])
        k_app = (n1415 ** 2) / (n1414 * n1515)
        print(temp, "    |", time, "     |", round(f1414, 3), "   |",
              round(f1415, 3), "   |", round(f1515, 3), "   |", round(k_app, 3))

if __name__ == "__main__":
     xef2_reaction()
     nitrogen_exchange()
     N14_N15_exchange()
