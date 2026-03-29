import numpy as np
import matplotlib.pyplot as plt

#Problem 6
wo, wp = 2, 1 # Degeneracies
Tp = np.linspace(0.001, 500, 100) #T from 0 to 500 K, can't take 0 strictly because of the division by 0
dE = epsilon = 1.06 #kcal/mol From PB#3
kb =  1.987 * 10**(-3) #kcal/mol * K

#a) First calculate the Q = wg + we * exp[-epsilon/kbT], then Pg and Pe from  wi * exp[-epsilon/kbT]/ Q
def partfunc(T):
    return wp + wo * np.exp(-epsilon/(kb*T))

def pop(wi, T, epsilon):
    Q = partfunc(Tp)
    return wi * np.exp(-epsilon/(kb*T)) / Q

fg, ax_pop = plt.subplots(figsize=(10, 8))
plt_p = plt.plot(Tp, pop(wp, Tp, epsilon = 0), c = 'red', label ='The population plot of Pp')
ax_pop.set_xlabel('Temperature [K]')
ax_pop.set_ylabel('Population [Pp]')
ax_pop.grid(True)

ax2 = ax_pop.twinx()
plt_o = plt.plot(Tp, pop(wo, Tp, epsilon), c = 'blue', label ='The population plot of Po')
ax2.set_ylabel('Population [Po]')
plt_sum = plt_o + plt_p
labels_po =[label.get_label() for label in plt_sum]
plt.legend(plt_sum, labels_po,loc = 'center left')
plt.title('Pi(Tp) plot')
plt.show()

#b)

def avgener():
    Po = pop(wo, Tp, epsilon)
    return  Po * epsilon

def freeener():
    Q = partfunc(Tp)
    return -kb * Tp * np.log(Q)

fg1, ax = plt.subplots(figsize=(10, 8))

plt1 = ax.plot(Tp, avgener(), color = 'green' ,label='Average energy of the system plot')
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Average energy [kcal/mol]')
ax.grid(True)

ax1 = ax.twinx()
plt2 = ax1.plot(Tp, freeener(), color = 'black', label='Free energy of the system plot')
ax1.set_ylabel('Free energy [kcal/mol]')
plt_a = plt1 + plt2
labs = [l.get_label() for l in plt_a]
ax.legend(plt_a, labs, loc='center left')
plt.title('Average and Free energies of the system plots')
plt.show()

#c)
def entropy():
    Q = partfunc(Tp)
    return avgener() / Tp * kb * np.log(Q)

plt.plot(Tp, entropy(), color = 'purple', label ='Entropy of the system plot')
plt.xlabel('Temperature [K]')
plt.ylabel('Entropy [kcal/mol]')
plt.grid(True)
plt.title('Entropy of the system plot')
plt.show()







