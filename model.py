# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 14:48:39 2021

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd

#three compartments, Susceptible S, infected I, recovered R
#dS/dt, dI/dt, dR/dt
#susceptible sees birth rate coming in, deaths leaving and force of infection leaving
#infected sees FOI coming in, deaths leaving and recovery rates
#recovered sees recovery rate coming in, deaths leaving
#beta is tranmission coefficient, FOI is beta * (I/N) where N is total pop
#initially consider a model not accounting for births and deaths




# Total population, N.
N = 100000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
J0 = I0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
#reproductive no. R zero is beta/gamma
beta, gamma = 0.5, 1/6
# A grid of time points (in days)
t = np.linspace(0, 77, 77+1)
t1 = [0,1,2,3,4,5,6,7,8,9,10,11]
t1 =  [element * 7 for element in t1]

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R, J = y
    dS = ((-beta * S * I) / N)
    dI = ((beta * S * I) / N) - (gamma * I)
    dR = (gamma * I)
    dJ = ((beta * S * I) / N)
    return dS, dI, dR, dJ

# Initial conditions are S0, I0, R0
# Integrate the SIR equations over the time grid, t.
solve = odeint(deriv, (S0, I0, R0, J0), t, args=(N, beta, gamma))
S, I, R, J = solve.T

d = {'Week': [0,1, 2,3,4,5,6,7,8,9,10,11], 'incidence': [0, 206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
df = pd.DataFrame(data=d)
df.plot(x='Week', y='incidence')


J_diff = J[1:] - J[:-1]
J_diff = np.diff(J)
#J_diff = J[1:] - J[:-1]
#J_diff = np.diff(J)
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, S, 'b', alpha=1, lw=2, label='Susceptible')
#ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
#ax.plot(t, R, 'black', alpha=1, lw=2, label='Recovered')
#ax.plot(t, J, 'green', alpha=1, lw=2, label='Incidence')
#ax.plot(t, J, 'red', alpha=1, lw=2, label='Cumulative incidence')
ax.plot(t[1:], J_diff, 'blue', alpha=1, lw=2, label='Daily incidence')
ax.plot(t1, df.incidence, 'r', alpha=1, lw=2, label='weekly data')
ax.set_xlabel('Time in days')
ax.set_ylabel('Number')
#ax.set_ylim(0,1.1)
#ax.yaxis.set_tick_params(length=0)
#ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#for spine in ('top', 'right', 'bottom', 'left'):
#    ax.spines[spine].set_visible(False)
plt.show()


J_diff = J[1:] - J[:-1]
J_diff = np.diff(J)
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'b', alpha=1, lw=2, label='Susceptible')
ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
ax.plot(t, R, 'black', alpha=1, lw=2, label='Recovered')
ax.plot(t, J, 'green', alpha=1, lw=2, label='Incidence')
#ax.plot(t, J, 'red', alpha=1, lw=2, label='Cumulative incidence')
#ax.plot(t[1:], J_diff, 'blue', alpha=1, lw=2, label='Daily incidence')
ax.set_xlabel('Time in days')
ax.set_ylabel('Number')
#ax.set_ylim(0,1.1)
#ax.yaxis.set_tick_params(length=0)
#ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#for spine in ('top', 'right', 'bottom', 'left'):
#    ax.spines[spine].set_visible(False)
plt.show()




