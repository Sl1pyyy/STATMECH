import numpy as np
import matplotlib.pyplot as plt

#Problems 1-2
Or = 85.3  # K
T = np.linspace(30, 300, 100) # K

list_T = []
list_f_para = []
list_f_ortho = []
list_cv_p = []
list_cv_o = []
list_cv_eq = []
list_cv_normal = []

for temp in T:
    q_p, e_p, e2_p = 0, 0, 0
    q_o, e_o, e2_o = 0, 0, 0

    for J in range(21):
        energy = J * (J + 1) * Or
        weight = (2 * J + 1) * np.exp(-energy / temp)

        if J % 2 == 0:  # Para-hydrogen (even J)
            q_p += 1 * weight
            e_p += 1 * energy * weight
            e2_p += 1 * (energy ** 2) * weight
        else:  # Ortho-hydrogen (odd J)
            q_o += 3 * weight
            e_o += 3 * energy * weight
            e2_o += 3 * (energy ** 2) * weight

    q_tot = q_p + q_o

    #Calculating heat capacities
    cv_p = ((e2_p / q_p) - (e_p / q_p) ** 2) / (temp ** 2) # J/(mol * K)
    cv_o = ((e2_o / q_o) - (e_o / q_o) ** 2) / (temp ** 2) # J/(mol * K)
    cv_eq = (((e2_p + e2_o) / q_tot) - ((e_p + e_o) / q_tot) ** 2) / (temp ** 2) # J/(mol * K)
    cv_normal = 0.25 * cv_p + 0.75 * cv_o # J/(mol * K)

    list_T.append(temp)
    list_f_para.append(q_p / q_tot)  # Fractions of para
    list_f_ortho.append(q_o / q_tot)  # Fractions of ortho
    list_cv_p.append(cv_p)
    list_cv_o.append(cv_o)
    list_cv_eq.append(cv_eq)
    list_cv_normal.append(cv_normal)

#Plotting 6 plots including experimental data
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

ax1.plot(list_T, list_f_para, label='Para-hydrogen fraction', color='blue')
ax1.plot(list_T, list_f_ortho, label='Ortho-hydrogen fraction', color='red')
ax1.set_title('Equilibrium fractions of para- and ortho-hydrogen')
ax1.set_ylabel('Mole fractions')
ax1.grid(True)
ax1.legend()

exp_T = [50, 100, 150, 200, 250, 298] # k
exp_Cv = [0.01, 0.15, 0.50, 0.78, 0.92, 0.98] # J/(mol * K)


ax2.plot(list_T, list_cv_o, label='Ortho-hydrogen heat capacity', color='red')
ax2.plot(list_T, list_cv_p, label='Para-hydrogen heat capacity', color='blue')
ax2.plot(list_T, list_cv_eq, label='Heat capacity of equilibrium mixture', color='green')
ax2.plot(list_T, list_cv_normal, label='Heat capacity of mixture (Frozen)', color='black', linestyle='--')
ax2.scatter(exp_T, exp_Cv, color='orange', label='Experimental heat capacity', zorder=5)

ax2.set_title('Heat capacity of hydrogen')
ax2.set_xlabel('Temperature [K]')
ax2.set_ylabel('Heat capacity [Cv/R]')
ax2.grid(True)
ax2.legend()
plt.tight_layout()
plt.show()