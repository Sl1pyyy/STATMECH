from math import factorial
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Problem 1

n = 7
nM = 3
nA = 2
nL = nS = 1

def probability1():
    """Function to calculate the probability of MAMMALS word to combine using the formulas:
     N = n!/ni and P(A) = N(A)/N in percents"""
    return 1/(factorial(n)/(factorial(nM) *factorial(nA) * factorial(nL)* factorial(nS))) * 100

problem1 = probability1()
if problem1:
    print(f'Probability that the word is “MAMMALS” : \n{problem1:2f}')

#Problem 2a

pdex2 = pd.read_excel('pbnum2.xlsx')
print(f'\nTable for problem 2: \n{pdex2}')

def norm_check():
    """Function to check whether the distribution is normalized or not"""
    return sum(pdex2['fRx'])

norm = norm_check()
if norm:
    print(f'\nThe value of column fRx: \n{norm}'
          f'\nColumn sum = 1, thus the  the distribution is normalized')

#Problem 2b
fig, ax1 = plt.subplots(figsize=(8, 4))

pdex2['cdf'] = np.cumsum(pdex2['fRx'])

ax1.bar(pdex2['Rx'], pdex2['fRx']/0.2, 0.1, label='Histogram of the f(Rg) and a plot of the F(Rg)')
ax1.grid(True)
ax1.set_xlabel('The radius of gyration,Rg')
ax1.set_ylabel('The probability distribution function,f(Rg)')
ax1.set_ylim(0, 1.2)

ax2 = ax1.twinx()
ax2.step(pdex2['Rx'],pdex2['cdf'],where='mid', color='black')
ax2.set_ylabel('The cumulative distribution function,F(Rg)')
plt.grid(True)
plt.show()

# Problem 2c
# The most probable value of Rg is Rg = 7.4 with probability = 0.205

rglist = []
frglist = []
for r in pdex2['Rx']:
    rglist.append(r)
for f in pdex2['fRx']:
    frglist.append(f)

rgavg = sum(rglisti * frglisti for rglisti,frglisti in zip(rglist,frglist))

print(f'The average value,Rgavg = {rgavg:2f}')

#Problem 2d
x = np.array(rglist)
p = np.array(frglist)

mean = np.sum(x * p)

#Formula sum( p * (x - mean)^2 )
variance = np.sum(p * (x - mean)**2)
std_dev = np.sqrt(variance)

#Formula sum( p * ((x - mean)/sigma)^3 )
skewness = np.sum(p * ((x - mean) / std_dev)**3)

#Formula sum( p * ((x - mean)/sigma)^4 )
kurtosis = np.sum(p * ((x - mean) / std_dev)**4)

print(f'\nVariance: {variance:.2f}')
print(f'\nSkewness: {skewness:.2f}')
print(f'\nKurtosis: {kurtosis:.2f}')

#Problem 3
##DATA for the exercise https://www.researchgate.net/publication/360116551_Pose_Classification_Using_Three-Dimensional_Atomic_Structure-Based_Neural_Networks_Applied_to_Ion_Channel-Ligand_Docking

N = 20
rmsd = np.array([3.53, 1.76, 3.95, 3.44, 3.83, 3.61, 3.42, 1.93, 2.59, 8.03, 4.29, 1.72, 4.11, 7.63, 4.45, 7.59, 5.31, 5.38, 3.61, 4.12])
bins_width = (np.max(rmsd) - np.min(rmsd))/len(rmsd)
print(f'\nBins width: {bins_width:.2f}')

bins = np.linspace(np.min(rmsd), np.max(rmsd), N + 1)
print(f'\n{bins}')

#ni calculations
counts, bin_edges = np.histogram(rmsd, bins=bins)

#Discrete probability (Pi = ni / N)
probabilities = counts / len(rmsd)

#Probability Density Function f(xi) = Pi / delta_x
bins_width = bin_edges[1] - bin_edges[0]
pdf_f_x = probabilities / bins_width

#Cumulative Distribution Function F(xi) = sum(Pi)
cdf_F_x = np.cumsum(probabilities)

dfex3 = pd.DataFrame({
    'bin_start': bin_edges[:-1],
    'bin_end': bin_edges[1:],
    'counts': counts,
    'PDF': pdf_f_x,
    'CDF': cdf_F_x
})

dfex3.to_excel('pbnum2ex3.xlsx')
print(f'\nExercise 3 functions:'
      f'\n{dfex3}')




