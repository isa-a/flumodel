# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 23:40:57 2021

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
#three compartments, Susceptible S, infected I, recovered R
#dS/dt, dI/dt, dR/dt
#susceptible sees birth rate coming in, deaths leaving and force of infection leaving
#infected sees FOI coming in, deaths leaving and recovery rates
#recovered sees recovery rate coming in, deaths leaving
#beta is tranmission coefficient, FOI is beta * (I/N) where N is total pop
#initially consider a model not accounting for births and deaths
# Total population, N.
N = 1000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 1, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
J0 = I0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 2/7, 1/7
# A grid of time points (in days)
t = np.linspace(0, 100,100)

  
empty = []
for i in range(100):
    empty.append(random.uniform(1.5, 2.5)*gamma)


# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R, J = y
    dS = ((-beta * S * I) / N)
    dI = ((beta * S * I) / N) - (gamma * I)
    dR = (gamma * I)
    dJ = ((beta * S * I) / N)
    return dS, dI, dR, dJ

solns =  []
for empt in empty:
    ces = odeint(deriv, (S0, I0, R0, J0), t, args=(N, empt, gamma))
    solns.append(ces)


J_diffs = []
for sol in solns:
    S, I, R, J = sol.T
    J_diffs.append(np.diff(J))
    

# plot all the solutions
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.set_xlabel('Time in days')
ax.set_ylabel('Number')
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
xcoords = [2.5, 97.5]
for J_diff in J_diffs:
    ax.plot(t[1:], J_diff, 'blue', alpha=1, lw=2)
    
def plot_percentile(data, percentile):
    percentiles = [np.percentile(val, percentile) for val in data]
    ax.plot(t, percentiles, 'blue', alpha=1, lw=2)

plot_percentile(solns, 97.5)
#plot_percentile(solns, 2.5)
#plot_percentile(solns, 0.5)

# plot without legend
plt.show()

# # Plot the data on three separate curves for S(t), I(t) and R(t)
# fig = plt.figure(facecolor='w')
# ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
# #ax.plot(t, S, 'b', alpha=1, lw=2, label='Susceptible')
# #ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
# #ax.plot(t, R, 'black', alpha=1, lw=2, label='Recovered')
# #ax.plot(t, J, 'green', alpha=1, lw=2, label='Incidence')
# #ax.plot(t, solns, 'green', alpha=1, lw=2, label='Incidence')
# ax.set_xlabel('Time in days')
# ax.set_ylabel('Number (1000s)')
# #ax.set_ylim(0,1.1)
# #ax.yaxis.set_tick_params(length=0)
# #ax.xaxis.set_tick_params(length=0)
# ax.grid(b=True, which='major', c='w', lw=2, ls='-')
# legend = ax.legend()
# legend.get_frame().set_alpha(0.5)
# #for spine in ('top', 'right', 'bottom', 'left'):
# #    ax.spines[spine].set_visible(False)
# plt.show()

for val in solns:
    ax.plot(t[1:], np.percentile(val,25), 'blue', alpha=1, lw=2)

for t in solns:
    ax.plot(t[1:], np.percentile(val,25), 'blue', alpha=1, lw=2)














