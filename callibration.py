# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 17:45:56 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from scipy.optimize import leastsq


def peak_infections(beta, days = 100):

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
    return (peak_infections(x) - 0.1) ** 2


res = minimize(residual, 0.5, method="Nelder-Mead")
print(res)

###############################################################################
##########                  WITH WEEKLY DATA
###############################################################################


t = np.arange(0,84,7)
d = {'Week': t, 'incidence': [0,206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
#d = {'Week': [time[7],time[14],time[21],time[28],time[35],time[42],time[49],time[56],time[63],time[70],time[77]], 'incidence': [206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
#d = {'Week': [1, 2,3,4,5,6,7,8,9,10,11], 'incidence': [206.1705794,2813.420201,11827.9453,30497.58655,10757.66954,7071.878779,3046.752723,1314.222882,765.9763902,201.3800578,109.8982006]}
df = pd.DataFrame(data=d)

def peak_infections(beta, df):

    # Weeks for which the ODE system will be solved
    weeks = df.Week.to_numpy()

    # Total population, N.
    N = 100000
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 10, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    J0 = I0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #reproductive no. R zero is beta/gamma
    gamma = 1/6 * 6 #rate should be in weeks now
    # A grid of time points 
    t = np.arange(0,84,7)

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

    return I/N

def residual(x, df):

    # Total population, N.
    N = 100000
    incidence = df.incidence.to_numpy()/N
    return np.sum((peak_infections(x, df) - incidence) ** 2)

x0 = 0.5
res = minimize(residual, x0, args=(df), method="Nelder-Mead").x
print(res)

best = leastsq(residual, x0,args=(df))
print(best)
