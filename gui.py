# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 17:26:08 2021

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import IntVar, DoubleVar
from PIL import ImageTk, Image


###############################################################################
def mainwindow():
    mainwindow = tk.Tk()
    mainwindow.geometry('350x350')
    
    tk.Label(mainwindow, text="enter parameters below").grid(row=1)
    
    getN = IntVar()
    geti0 = IntVar()
    getr0 = IntVar()
    #getbeta = IntVar()
    #getgamma = IntVar()
    
    getbeta = DoubleVar()
    getgamma = DoubleVar()

    
    tk.Label(mainwindow, text="N").grid(row=2)
    tk.Label(mainwindow, text="i0").grid(row=3)
    tk.Label(mainwindow, text="r0").grid(row=4)
    tk.Label(mainwindow, text="beta").grid(row=5)
    tk.Label(mainwindow, text="gamma").grid(row=6)
    
    e1 = tk.Entry(mainwindow,textvariable = getN).grid(row=2, column=1)
    e2 = tk.Entry(mainwindow,textvariable = geti0).grid(row=3, column=1)
    e3 = tk.Entry(mainwindow,textvariable = getr0).grid(row=4, column=1)
    e4 = tk.Entry(mainwindow,textvariable = getbeta).grid(row=5, column=1)
    e5 = tk.Entry(mainwindow,textvariable = getgamma).grid(row=6, column=1)
    
    solve = tk.Button(mainwindow, text='solve!', command=lambda: [values()]).grid(row=7, column=1, sticky=tk.W, pady=4)
    
    
    
    def values():

        global readN
        global readi0
        global readr0
        global readbeta
        global readgamma
        
        readN = getN.get()
        readi0 = geti0.get()
        readr0 = getr0.get()
        readbeta = (getbeta.get())
        readgamma =(getgamma.get())
        
        
        #print('ppppppppppppp', readbeta,readgamma)
        
        intN = int(readN)
        inti0 = int(readi0)
        intr0 = int(readr0)
        intbeta = float(readbeta)       #ENTER AS DECIMAL
        intgamma = float(readgamma)     #ENTER AS DECIMAL
        
        
        # Total population, N.
        N = readN
        # Initial number of infected and recovered individuals, I0 and R0.
        I0, R0 = readi0, readr0
        # Everyone else, S0, is susceptible to infection initially.
        S0 = readN - readi0 - readr0
        # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
        beta, gamma = readbeta, readgamma
        # A grid of time points (in days)
        t = np.linspace(0, 160, 160)

        # The SIR model differential equations.
        def deriv(y, t, N, I0, R0):
            S, I, R = y
            dS = ((-beta * S * I) / N)
            dI = ((beta * S * I) / N) - (gamma * I)
            dR = (gamma * I)
            return dS, dI, dR
        
        # Initial conditions are S0, I0, R0
        # Integrate the SIR equations over the time grid, t.
        solve = odeint(deriv, (S0, I0, R0), t, args=(N, beta, gamma))
        S, I, R = solve.T
        
        # Plot the data on three separate curves for S(t), I(t) and R(t)
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
        ax.plot(t, S/1000, 'blue', alpha=1, lw=2, label='Susceptible')
        ax.plot(t, I/1000, 'r', alpha=1, lw=2, label='Infected')
        ax.plot(t, R/1000, 'black', alpha=1, lw=2, label='Recovered')
        ax.set_xlabel('Time in days')
        ax.set_ylabel('Number (1000s)')
        #ax.set_ylim(0,1.2)
        ax.yaxis.set_tick_params(length=0)
        ax.xaxis.set_tick_params(length=0)
        ax.grid(b=True, which='major', c='w', lw=2, ls='-')
        legend = ax.legend()
        legend.get_frame().set_alpha(0.5)
        for spine in ('top', 'right', 'bottom', 'left'):
            ax.spines[spine].set_visible(False)
        plt.title('Influenza dynamics')
        plt.savefig('N=' + str(readN) + ', i0=' + str(readi0) + 
                    ', r0=' + str(readr0) + ', beta=' + str(readbeta) + ', gamma=' + str(readgamma) + '.png')
        #plt.show()
        
        def plotwindow():
            #This creates the main window of an application
            window = tk.Toplevel()
            window.title("Join")
            window.geometry("500x400")
            window.configure(background='grey')
            path = ('N=' + str(readN) + ', i0=' + str(readi0) + 
                                ', r0=' + str(readr0) + ', beta=' + str(readbeta) + ', gamma=' + str(readgamma) + '.png')
            
            #Creates a Tkinter-compatible photo image, which can be used everywhere Tkinter expects an image object.
            img = ImageTk.PhotoImage(Image.open(path))
            
            #The Label widget is a standard Tkinter widget used to display a text or image on the screen.
            panel = tk.Label(window, image = img)
            
            #The Pack geometry manager packs widgets in rows or columns.
            panel.pack(side = "bottom", fill = "both", expand = "yes")
            
            #Start the GUI
            window.mainloop()
        plotwindow()


    mainwindow.mainloop()
    
mainwindow()


