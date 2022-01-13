# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 17:45:56 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from  scipy.optimize import root
import pandas as pd


def peak_infections(beta, days = 100):

    # Total population, N.
    N = 100000
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 10, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    J0 = I0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
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


betas = np.linspace(0,1,101,endpoint = True)
peak_inf = [peak_infections(b) for b in betas]
plt.plot(betas, peak_inf)
plt.plot(betas, 0.1*np.ones(len(betas)))



root(lambda b: peak_infections(b)-0.3049758655, x0 = 0.5).x


data = pd.read_csv('data.csv')
x = data['Week']
y = data['incidence']
plt.plot(x,y)


