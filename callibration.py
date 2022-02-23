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
d = {'Week': [t[0], t[7],t[14],t[21],t[28],t[35],t[42],t[49],t[56],t[63],t[70],t[77]], 
     'incidence': [0, 206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,
                   7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
df = pd.DataFrame(data=d)
#d = {'Week': t, 'incidence': [0,206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
#df = pd.DataFrame(data=d)

def peak_infections(beta, df):
 
    # Weeks for which the ODE system will be solved
    #weeks = df.Week.to_numpy()

    # Total population, N.
    N = 100000
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 10, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    J0 = I0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #reproductive no. R zero is beta/gamma
    gamma = 1/6 #rate should be in weeks now
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

x0 = 0.5
res = minimize(residual, x0, args=(df), method="Nelder-Mead", options={'fatol':1e-04}).x
print(res)

plt.plot(d['Week'], df.incidence.to_numpy()/100000, label="Real data")
plt.plot(d['Week'], peak_infections(.53, df), label="Model with 0.72")
#plt.plot(d['Week'], peak_infections(.42, df), label="Model with .42")
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

