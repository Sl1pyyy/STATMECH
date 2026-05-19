import matplotlib.pyplot as plt
import math

R = 8.314 * 10 ** (-3)  # kJ / (mol * K)
H_cm = 1.4388 # K

#Problem 1

def plot_xef2_reaction_energy():
    R = 8.314  # J/(mol*K)
    h_cm = 1.4388  # K

    #Electronic energy change in kJ/mol
    # -267.8 - (-154.8) = -113.0 kJ/mol
    E_elec = -113.0

    wavenumbers_xef2 = [214, 214, 515, 555]  # 4 total modes (214 is degenerate)
    wavenumbers_f2 = [894]  # 1 mode

    theta_xef2 = []
    for wn in wavenumbers_xef2:
        theta_xef2.append(wn * h_cm)

    theta_f2 = []
    for wn in wavenumbers_f2:
        theta_f2.append(wn * h_cm)

    temperature = []
    U = []

    for T in range(300, 2010, 10):

        u_trans = -1.5 * R * T

        u_vib_xef2 = 0
        for theta in theta_xef2:
            u_vib_xef2 = u_vib_xef2 + (R * theta) / (math.exp(theta / T) - 1)

        u_vib_f2 = 0
        for theta in theta_f2:
            u_vib_f2 = u_vib_f2 + (R * theta) / (math.exp(theta / T) - 1)

        U_th = u_trans + u_vib_xef2 - u_vib_f2
        E_total = E_elec + (U_th / 1000.0)
        temperature.append(T)
        U.append(E_total)

    plt.figure(figsize=(8, 5))
    plt.plot(temperature, U, color="skyblue", linewidth=2.5)
    plt.title("Total reaction energy change (Ureact.)", fontsize=12)
    plt.xlabel("Temperature [K]", fontsize=11)
    plt.ylabel("Ureact. [kJ/mol]", fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.show()

#Problem 2
def plot_nitrogen_isotope_exchange():
    m14 = 14.003074
    m15 = 15.000108

    #Given characteristic vibrational temperature for 14N2
    theta_14_14 = 3374.0  # K

    mu_14_14 = m14 / 2.0
    mu_15_15 = m15 / 2.0
    mu_14_15 = (m14 * m15) / (m14 + m15)

    theta_15_15 = theta_14_14 * math.sqrt(mu_14_14 / mu_15_15)
    theta_14_15 = theta_14_14 * math.sqrt(mu_14_14 / mu_14_15)

    dZPE = (2.0 * theta_14_15 - theta_14_14 - theta_15_15) / 2.0

    mass_factor = 2.0 * (m14 + m15) / math.sqrt(m14 * m15)

    t_values = list(range(200, 5200, 200))

    log_K_list = []
    x_labels = []
    x_positions = []

    for idx, T in enumerate(t_values):
        zpe_term = math.exp(-dZPE / T)

        num = (1.0 - math.exp(-theta_14_14 / T)) * (1.0 - math.exp(-theta_15_15 / T))
        den = (1.0 - math.exp(-theta_14_15 / T)) ** 2
        vib_term = num / den

        K = mass_factor * zpe_term * vib_term

        log_K_list.append(math.log10(K))

        x_positions.append(idx)
        inv_T = 1000.0 / T
        if inv_T == int(inv_T):
            x_labels.append(f"{int(inv_T)}")
        else:
            x_labels.append(f"{inv_T:.3g}")

    plt.figure(figsize=(11, 5))
    plt.plot(x_positions, log_K_list, color="green", linewidth=2, label=" Arrhenius plot log K = f(1/T)")
    plt.axhline(y=math.log10(4.0), color="coral", linestyle="--", label="lg(4) = 0.602 asymptotic high-temperature thermodynamic limit")
    plt.title("Arrhenius plot log K = f(1/T)", fontsize=12)
    plt.xlabel("$1000/T\ (\mathrm{K}^{-1})$", fontsize=11)
    plt.ylabel("log_(10) K", fontsize=11)
    plt.xlim(-0.5, len(x_positions) - 0.5)
    plt.ylim(0.6012, 0.6021)
    plt.grid(True, linestyle="-", alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.show()


#Problem 3

def isotope_exchange():
    #Raw experimental data points: (Temperature, time, n_14_15, n_15_15)
    data_points = [
        # 465 °C
        {"T_C": 465, "time": 2, "n_14_15": 100.5, "n_15_15": 51.50},
        {"T_C": 465, "time": 17, "n_14_15": 125.8, "n_15_15": 39.00},
        {"T_C": 465, "time": 37, "n_14_15": 163.0, "n_15_15": 22.60},
        {"T_C": 465, "time": 60, "n_14_15": 194.0, "n_15_15": 14.00},
        # 500 °C
        {"T_C": 500, "time": 2, "n_14_15": 108.0, "n_15_15": 48.00},
        {"T_C": 500, "time": 8, "n_14_15": 142.0, "n_15_15": 38.00},
        {"T_C": 500, "time": 21, "n_14_15": 184.0, "n_15_15": 20.00},
        {"T_C": 500, "time": 26, "n_14_15": 197.5, "n_15_15": 17.40},
        {"T_C": 500, "time": 35, "n_14_15": 202.0, "n_15_15": 13.40},
        {"T_C": 500, "time": 46, "n_14_15": 207.0, "n_15_15": 13.25}
    ]

    current_temp = None

    for pt in data_points:
        if pt["T_C"] != current_temp:
            current_temp = pt["T_C"]
            print(f"\n--- Fractions and K — {current_temp} °C ---")
            print(
                f"{'t (h)':<6}{'n(14N2)':<10}{'n(15N2)':<10}{'n(14N15N)':<12}{'f(14N2)%':<11}{'f(15N2)%':<11}{'f(14N15N)%':<13}{'K_app':<6}")
            print("-" * 85)

        n_14_14 = 1000.0 - (pt["n_14_15"] / 2.0)
        n_total = n_14_14 + pt["n_15_15"] + pt["n_14_15"]
        f_14_14 = (n_14_14 / n_total) * 100.0
        f_15_15 = (pt["n_15_15"] / n_total) * 100.0
        f_14_15 = (pt["n_14_15"] / n_total) * 100.0
        K_app = (f_14_15 / 100.0) ** 2 / ((f_14_14 / 100.0) * (f_15_15 / 100.0))
        print(
            f"{pt['time']:<6}{n_14_14:<10.1f}{pt['n_15_15']:<10.2f}{pt['n_14_15']:<12.1f}{f_14_14:<11.2f}{f_15_15:<11.2f}{f_14_15:<13.2f}{K_app:<6.3f}")

if __name__ == "__main__":
     plot_xef2_reaction_energy()
     plot_nitrogen_isotope_exchange()
     isotope_exchange()
