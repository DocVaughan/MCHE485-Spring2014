#! /usr/bin/env python

##########################################################################################
# mass_spring_damper_InitialCond.py
#
# Script to analyze a spring-mass-damper system response to initial conditions
#
# Uses the control systems module that can be found here:
#   https://www.cds.caltech.edu/~murray/wiki/Control_Systems_Library_for_Python
#
# Created: 01/23/14
#   - Joshua Vaughan
#   - joshua.vaughan@louisiana.edu
#   - http://www.ucs.louisiana.edu/~jev9637
#
# Modified:
#   *
#
##########################################################################################

# Simple mass-spring system with Force input
#                   +---> X
#                   |
#   /|     k     +-----+
#   /|---/\/\/---|     |
#   /|           |  M  +=====> F
#   /|-----]-----|     |
#   /|     c     +-----+

import numpy as np              # Grab all of the NumPy functions with nickname np
from matplotlib.pyplot import * # Grab MATLAB plotting functions
import control                  # import the control system functions

# Uncomment to use LaTeX to process the text in figure
# from matplotlib import rc
# rc('text',usetex=True)

m = 1.0                 # kg
k = (2.*np.pi)**2.      # N/m (Selected to give an undamped wn of 1Hz)
wn = np.sqrt(k/m)       # Natural Frequency (rad/s)

z = 0.1                 # Define a desired damping ratio
c = 2*z*wn*m            # calculate the damping coeff. to create it (N/(m/s))


# Define the system to use in simulation - in transfer function form here
num = [1/m]
den = [1,2.*z*wn,wn**2]

sys = control.tf(num,den)


# Set up simulation parameters
t = np.linspace(0,5,500)            # time for simulation, 0-5s with 500 points in-between

F = np.zeros_like(t)

# Define the initial conditions x_dot(0) = 0, x(0) = 1
x0 = [0.,1.]

# run the simulation - utilize the built-in initial condition response function
[T,yout] = control.initial_response(sys,t,x0)

# Make the figure pretty, then plot the results
#   "pretty" parameters selected based on pdf output, not screen output
#   Many of these setting could also be made default by the .matplotlibrc file
fig = figure(figsize=(6,4))
ax = gca()
subplots_adjust(bottom=0.17,left=0.17,top=0.96,right=0.96)
setp(ax.get_ymajorticklabels(),family='CMU Serif',fontsize=18)
setp(ax.get_xmajorticklabels(),family='CMU Serif',fontsize=18)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.grid(True,linestyle=':',color='0.75')
ax.set_axisbelow(True)

xlabel('Time (s)',family=   ',fontsize=22,weight='bold',labelpad=5)
ylabel('Position (m)',family='CMU Serif',fontsize=22,weight='bold',labelpad=10)
ylim(-1.,1.)

plot(T,yout,color="blue",linewidth=2)

# Save the figure as a high-res pdf in the current folder
if z < 1.:
    savefig('Underdamped_Resp.pdf')
else:
    savefig('Overdamped_Resp.pdf')

# show the figure
show()


# We could also use the solution to the differential equation 
#   that we developed in class to plot the response. 

# Define x(t) for the underdamped case
x = np.exp(-z*wn*t)*(x0[1]*np.cos(wd*t) + (z*wn*x0[1] + x0[0])/wd * np.sin(wd*t))

# Now just plot the response

# Make the figure pretty, then plot the results
#   "pretty" parameters selected based on pdf output, not screen output
#   Many of these setting could also be made default by the .matplotlibrc file
fig = figure(figsize=(6,4))
ax = gca()
subplots_adjust(bottom=0.17,left=0.17,top=0.96,right=0.96)
setp(ax.get_ymajorticklabels(),family='serif',fontsize=18)
setp(ax.get_xmajorticklabels(),family='serif',fontsize=18)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.grid(True,linestyle=':',color='0.75')
ax.set_axisbelow(True)

xlabel('Time (s)',family='serif',fontsize=22,weight='bold',labelpad=5)
ylabel('Position (m)',family='serif',fontsize=22,weight='bold',labelpad=10)
ylim(-1.,1.)

plot(t,x,color="blue",linewidth=2)
