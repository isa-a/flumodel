# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:44:23 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd

#seven compartments, Susceptible S, vaccinated s1a, s1b, double dosed s2a, s2b, infected I, recovered R
#dS/dt, dS1a/dt, dS1b/dt, dS2a/dt, dS2b/dt, dI/dt, dR/dt
#susceptible sees birth rate coming in, deaths leaving and force of infection leaving
#infected sees FOI coming in, deaths leaving and recovery rates
#recovered sees recovery rate coming in, deaths leaving
#beta is tranmission coefficient, FOI is beta * (I/N) where N is total pop
#initially consider a model not accounting for births and deaths

# Total population, N.
N = 1000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
#first vaccinated, s1a
S1a = 0.3*(S0)
#second, s1b
S1b = 0.4*(S1a)
#s2a
S2a = 0.1*(S1a)
#s2b
S2b = 0.2*(S1b)
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 1/3, 1/7
x1a, x1b, x2a, x2b = 0.1, 0.13, 0.15, 0.21
v, va, vb = 0.3, 0.2, 0.23
# A grid of time points (in days)
t = np.linspace(0, 100, 100+1)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma, x1a, x1b, x2a, x2b, v, va, vb):
    S, S1a, S1b, S2a, S2b, I, R = y
    
    dS = ((-beta * S * I) / N)
    
    dS1a = v*S - (((beta * S1a * I) / N) * x1a) - t*S1a - va*S1a
    dS1b = t*S1a - (((beta * S1b * I) / N) * x1b) - vb*S1b
    
    dS2a = va*S1a - (((beta * S2a * I) / N) * x2a) - t*S2a
    dS2b = t*S2a - (((beta * S2b * I) / N) * x2b) + vb*S1b
    
    dI = ((beta * S * I) / N) + (((beta * S1a * I) / N) * x1a) + (((beta * S2a * I) / N) * x2a) + (((beta * S2b * I) / N) * x2b) + (((beta * S1b * I) / N) * x1b) - (gamma * I)
    
    dR = (gamma * I)

    return dS, dS1a, dS1b, dS2a, dS2b, dI, dR

# Initial conditions are S0, I0, R0
# Integrate the SIR equations over the time grid, t.
solve = odeint(deriv, (S0, I0, R0, S1a, S1b, S2a, S2b), t, args=(N, beta, gamma, x1a, x1b, x2a, x2b, v, va, vb))
S, S1a, S1b, S2a, S2b, I, R = solve.T

#J_diff = J[1:] - J[:-1]
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, S/1000, 'b', alpha=1, lw=2, label='Susceptible')
ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
#ax.plot(t, R/1000, 'black', alpha=1, lw=2, label='Recovered')
#ax.plot(t, J/1000, 'green', alpha=1, lw=2, label='Incidence')
#ax.plot(t, J, 'red', alpha=1, lw=2, label='Cumulative incidence')
#ax.plot(t[1:], J_diff, 'blue', alpha=1, lw=2, label='Daily incidence')
ax.set_xlabel('Time in days')
ax.set_ylabel('Number (1000s)')
#ax.set_ylim(0,1.1)
#ax.yaxis.set_tick_params(length=0)
#ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#for spine in ('top', 'right', 'bottom', 'left'):
#    ax.spines[spine].set_visible(False)
plt.show()