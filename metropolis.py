# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 02:09:52 2022

@author: ISA
"""

import numpy as np
import random
import matplotlib.pyplot as plt


def normaldist(x,mew,sigma): #'normal distribution formula
    num = np.exp((-(x-mew)**2)/(2*sigma**2))
    den = sigma * np.sqrt(2*np.pi)
    return num/den

def uniformdist(p): #uniform distribution
    unif = random.uniform(0,1)
    if unif>=p:
        return False
    else:
        return True
    
def gaussian_mcmc(jump,mew,sigma):
    states = []
    burn_in = int(jump*0.2)
    current = random.uniform(-5*sigma+mew,5*sigma+mew)
    for i in range(jump):
        states.append(current)
        movement = random.uniform(-5*sigma+mew,5*sigma+mew)
        
        curr_prob = normaldist(x=current,mew=mew,sigma=sigma)
        move_prob = normaldist(x=movement,mew=mew,sigma=sigma)
        
        acceptance = min(move_prob/curr_prob,1)
        if uniformdist(acceptance):
            current = movement
    return states[burn_in:]
    
lines = np.linspace(-3,3,1000)
normal_curve = [normaldist(l,mew=0,sigma=1) for l in lines]
dist = gaussian_mcmc(100_000,mew=0,sigma=1)
plt.hist(dist,density=1,bins=20) 
plt.plot(lines,normal_curve)