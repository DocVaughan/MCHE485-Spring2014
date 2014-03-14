#! /usr/bin/env python 

# Script to generate compare viscous damping to friction forec

import numpy as np
from matplotlib.pyplot import * 


x_dot = np.linspace(-1,1,500)
x_dot_fr = np.linspace(-1,1,500)

# viscous damping force
c = 0.1
F_vis = c*x_dot

# friction
u = 0.1
N = 1.
F_fr = u*N*np.sign(x_dot)

# mask the discontinuity
pos = np.where(np.abs(x_dot_fr) <= 1e-2)
F_fr[pos] = np.nan
x_dot_fr[pos] = np.nan

#-----  Copy from here down into your code, replacing items as needed ----------------
#
# Set the plot size - 3x2 aspect ratio is best
fig = figure(figsize=(6,4))
ax = gca()
subplots_adjust(bottom=0.17,left=0.17,top=0.96,right=0.96)

# Change the axis units to CMU Serif
setp(ax.get_ymajorticklabels(),family='CMU Serif',fontsize=18)
setp(ax.get_xmajorticklabels(),family='CMU Serif',fontsize=18)

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Turn on the plot grid and set appropriate linestyle and color
ax.grid(True,linestyle=':',color='0.75')
ax.set_axisbelow(True)

# Define the X and Y axis labels
xlabel('Relative Velocity, $\dot{x}$',family='CMU Serif',fontsize=22,weight='bold',labelpad=5)
ylabel('Force (N)',family='CMU Serif',fontsize=22,weight='bold',labelpad=10)

plot(x_dot,F_vis,linewidth=2,linestyle="--",label=r'Viscous')
plot(x_dot_fr,F_fr,linewidth=2,label=r'Coulomb')

# uncomment below and set limits if needed
# xlim(0,5)
ylim(-.11,.11)
xticks([0])
yticks([-0.1, 0, 0.1],['$-\mu N$', 0, '$\mu N$'])

# Create the legend, then fix the fontsize
leg = legend(loc='lower right', fancybox=True, borderaxespad = 0.5)
ltext  = leg.get_texts()
setp(ltext,family='CMU Serif',fontsize=16)

# Adjust the page layout filling the page using the new tight_layout command
tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
savefig('FrictionViscous_comparison.pdf',dpi=300) 
savefig('FrictionViscous_comparison.pdf',dpi=300) 

show()