# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 17:45:56 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
import pandas as pd
from scipy.optimize import leastsq

###############################################################################
##########                  WITH WEEKLY DATA
###############################################################################


#t = np.arange(0,84,7)
t = np.linspace(0, 77, 77+1)
d = {'Week': [t[0],t[7],t[14],t[21],t[28],t[35],t[42],t[49],t[56],t[63],t[70],t[77]], 
     'incidence': [10,206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,
                   7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
df = pd.DataFrame(data=d)
#d = {'Week': t, 'incidence': [0,206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
#df = pd.DataFrame(data=d)

def peak_infections(x, df):
 
    # Weeks for which the ODE system will be solved
    #weeks = df.Week.to_numpy()

    # Total population, N.
    N = 100000
    # Initial number of infected and recovered individuals, I0 and R0.
    #I0, R0 = 10, 0
    # Everyone else, S0, is susceptible to infection initially.
    R0 = 0
    beta = x[0]
    I0 = x[1]
    gamma = x[2]
    S0 = N - I0 - R0
    J0 = I0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #reproductive no. R zero is beta/gamma
    #gamma = 1/6 #rate should be in weeks now
    # A grid of time points 
    times = np.arange(0,84,7)

    # The SIR model differential equations.
    def deriv(y, times, N, beta, gamma):
        S, I, R, J = y
        dS = ((-beta * S * I) / N)
        dI = ((beta * S * I) / N) - (gamma * I)
        dR = (gamma * I)
        dJ = ((beta * S * I) / N) #incidence
        return dS, dI, dR, dJ

    # Initial conditions are S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    solve = odeint(deriv, (S0, I0, R0, J0), times, args=(N, beta, gamma))
    S, I, R, J = solve.T

    return I/N

def residual(x, df):

    # Total population,  N.
    StartingPop = 100000
    incidence = df.incidence.to_numpy()/StartingPop
    return np.sum((peak_infections(x,df) - incidence) ** 2)


x0 = [0.5, 10, 1/6] #beta, i0, gamma
res = minimize(residual, x0, args=(df), method="Nelder-Mead", options={'fatol':1e-04}).x
print(res)

initial_infecteds = np.arange(1,31,1)
initial_beta = np.arange(0,1,1/10)
#initial_gamma = np.arange(1/10,1/6,1/30)
       
def fit_infecteds(i):
    guesses = [0.5, i, 1/6] #beta, i0, gamma
    res = minimize(residual, guesses, args=(df), method="Nelder-Mead", options={'fatol':1e-04}).x
    print(res, f'for a value of I0={i}')

for i in initial_infecteds:
    fit_infecteds(i)

def fit_beta(j):
    guesses = [j, 10, 1/6] #beta, i0, gamma
    res = minimize(residual, guesses, args=(df), method="Nelder-Mead", options={'fatol':1e-04}).x
    print(res, f'for a value of {j}')

for j in initial_beta:
    fit_beta(j)


fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, S, 'b', alpha=1, lw=2, label='Susceptible')
#ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
#ax.plot(t, R, 'black', alpha=1, lw=2, label='Recovered')
#ax.plot(t, J, 'green', alpha=1, lw=2, label='Incidence')
#ax.plot(t, J, 'red', alpha=1, lw=2, label='Cumulative incidence')
ax.plot(d['Week'], df.incidence.to_numpy(), label="Real data")
ax.plot(d['Week'], peak_infections([0.66,1.41,0.21], df)*100000, label="Model with optimal values")
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


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from scipy.optimize import leastsq

#t = np.arange(0,84,7)
t = np.linspace(0, 77, 77+1)
d = {'Week': [t[7],t[14],t[21],t[28],t[35],t[42],t[49],t[56],t[63],t[70],t[77]], 'incidence': [206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
df = pd.DataFrame(data=d)
#d = {'Week': t, 'incidence': [0,206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
#df = pd.DataFrame(data=d)

def peak_infections(beta, df):
 
    # Weeks for which the ODE system will be solved
    #weeks = df.Week.to_numpy()

    # Total population, N.
    N = 100000
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 1.41351201, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    J0 = I0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #reproductive no. R zero is beta/gamma
    gamma = 0.21620529 #rate should be in weeks now
    # A grid of time points 
    t7 = np.arange(7,84,7)

    # The SIR model differential equations.
    def deriv(y, t7, N, beta, gamma):
        S, I, R, J = y
        dS = ((-beta * S * I) / N)
        dI = ((beta * S * I) / N) - (gamma * I)
        dR = (gamma * I)
        dJ = ((beta * S * I) / N)
        return dS, dI, dR, dJ

    # Initial conditions are S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    solve = odeint(deriv, (S0, I0, R0, J0), t7, args=(N, beta, gamma))
    S, I, R, J = solve.T

    return np.max(I)/N

def residual(x, df):

    # Total population,  N.
    N = 100000
    incidence = df.incidence.to_numpy()/N
    return np.sum((peak_infections(x, df) - incidence) ** 2)

x0 = 0.66411585
res = minimize(residual, x0, args=(df), method="Nelder-Mead").x
print(res)


fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, S, 'b', alpha=1, lw=2, label='Susceptible')
#ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
#ax.plot(t, R, 'black', alpha=1, lw=2, label='Recovered')
#ax.plot(t, J, 'green', alpha=1, lw=2, label='Incidence')
#ax.plot(t, J, 'red', alpha=1, lw=2, label='Cumulative incidence')
ax.plot(d['Week'], df.incidence.to_numpy()/100000, label="Real data")
ax.plot(d['Week'], peak_infections(0.53, df)*100000, label="Model with 0.66")
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

import matplotlib.pyplot as plt
plt.plot(d['Week'], df.incidence.to_numpy()/100000, label="Real data")
plt.plot(d['Week'], peak_infections(.72, df), label="Model with 0.72")
plt.plot(d['Week'], peak_infections(.42, df), label="Model with .42")
plt.legend()

# =============================================================================
# res2 = leastsq(residual, x0,args=(df))
# print(res2)
# 
# results = minimize_scalar(residual,(0.4, 0.5),args=(df))
# print(results)
# results['fun']
# 
# =============================================================================

###############################################################################
##########                  --------------------------
###############################################################################


def peak_infections_days(beta, days = 100):

    # Total population, N.
    N = 1000
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 10, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    J0 = I0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #reproductive no. R zero is beta/gamma
    gamma = 1/7
    # A grid of time points (in days)
    t = np.linspace(0, days, days + 1)

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

    return np.max(I)/N


#scipy.optimize

data = pd.read_csv('data.csv')
x = data['Week']
y = data['incidence']
plt.plot(x,y)



def residual(x):
    return (peak_infections_days(x) - 0.1) ** 2


res = minimize(residual, 0.5, method="Nelder-Mead").x
print(res)

best = leastsq(residual, x0)
print(best)

