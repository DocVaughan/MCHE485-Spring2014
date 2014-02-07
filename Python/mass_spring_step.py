#! /usr/bin/env python 

##########################################################################################
# mass_spring_step.py
#
# Script to analyze a spring-mass-damper system response to step input in position
# 
# Created: 2/4/13 
#   - Joshua Vaughan 
#   - joshua.vaughan@louisiana.edu
#   - http://www.ucs.louisiana.edu/~jev9637
#
# Modified:
#   *
#
##########################################################################################

# The baseline system, a simple mass-spring-damper system
#
#    +---> Y       +---> X
#    |             |
#    |    k     +-----+
#    +---/\/\---|     |
#    |          |  M  |
#    +----]-----|     |
#         c     +-----+
#
# Step input in Y, response of the mass is X

from numpy import *             # Grab all of the NumPy functions
from matplotlib.pyplot import * # Grab MATLAB plotting functions
import control                  # import the control system functions

# Uncomment to use LaTeX to process the text in figure
# from matplotlib import rc
# rc('text',usetex=True)

m = 1.;             # kg
k = (2.*pi)**2;     # N/m (Selected to give an undamped wn of 1Hz)
wn = sqrt(k/m)      # Natural Frequency (rad/s)

z = 0.25;           # Define a desired damping ratio
c = 2*z*wn*m;       # calculate the damping coeff. to create it (N/(m/s))


# Define the system to use in simulation - in transfer function form here
num = [2.*z*wn,wn**2];
den = [1,2.*z*wn,wn**2];

sys = control.tf(num,den);
       

# Set up simulation parameters
t = np.linspace(0,3,500)            # Time for simulation, 0-5s with 500 points in-between



# run the simulation - utilize the built-in initial condition response function
[T,yout] = control.step_response(sys,t)

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

xlabel('Time (s)',family='CMU Serif',fontsize=22,weight='bold',labelpad=5)
ylabel('Position (m)',family='CMU Serif',fontsize=22,weight='bold',labelpad=10)

plot(T,yout,color="blue",linewidth=2)

# save the figure as a high-res pdf in the current folder
savefig('Step_Resp.pdf',dpi=300)

# show the figure
show()


#-----  We could also manually input the step response -----
# Let's start it at t=0.5s for
U = zeros(500)  # Define an array of all zeros
U[50:] = 1      # Make all elements of this array index>50 = 1 (all after 0.5s)

# run the simulation - utilize the built-in initial condition response function
[T_man,yout_man,xout_man] = control.forced_response(sys,t,U)

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

xlabel('Time (s)',family='CMU Serif',fontsize=22,weight='bold',labelpad=5)
ylabel('Position (m)',family='CMU Serif',fontsize=22,weight='bold',labelpad=10)

plot(T_man,U,color='red',linewidth=2,linestyle='--',label='Step Input')
plot(T_man,yout_man,color='blue',linewidth=2,label='Response')

leg = legend(loc='upper right', fancybox=True)
ltext  = leg.get_texts() 
setp(ltext,family='CMU Serif',fontsize=16)

# save the figure as a high-res pdf in the current folder
savefig('Step_Resp_Manual.pdf',dpi=600)

# show the figure
show()



