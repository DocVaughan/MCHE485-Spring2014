#! /usr/bin/env python 

##########################################################################################
# Script to animate a numerical approximation of a Fourier Expansion
#
#                         ----- THIS IS NOT AN FFT!!! -----
# 
# Created: 3/125/13 
#   - Joshua Vaughan 
#   - joshua.vaughan@louisiana.edu
#   - http://www.ucs.louisiana.edu/~jev9637
#
# Modified:
#   * 
#
##########################################################################################

import numpy as np                      # Grab all of the NumPy functions
from matplotlib import pyplot as plt    # Grab MATLAB plotting functions
import control                          # import the control system functions
import matplotlib.animation as animation


#--------- Input your function to examine here --------
#
# This is the square wave from the Fall 2012 midterm example we did in class
# t = np.linspace(0,4,4000)       # define the time to look at, easiest to just choose 1 period
# w0 = 2*np.pi/t[-1]              # define the fundamental frequency (here, I know t(end)=tau)
# tau_0 = 2*np.pi/w0              # define fundamental period based on w0
# #
# y = (t>-4) - 2*(t>-2) + 2*(t>0) - 2*(t>2) + 2*(t>4) - 2*(t>6) + 2*(t>8) - (t>10)

# This is the "trapezoid" wave example we worked in class
t = np.linspace(0,4,4000)       # define the time to look at, easiest to just choose 1 period
w0 = 2*np.pi/t[-1]                 # define the fundamental frequency (here, I know t(end)=tau)
tau_0 = 2*np.pi/w0                 # define fundamental period based on w0
#
F0 = 1
y = np.zeros((len(t),))

for ii in range(len(t)):
  if t[ii] <= tau_0/3:
      y[ii] = 3*F0/tau_0*t[ii]
  elif t[ii] <= 2*tau_0/3:
      y[ii] = F0
  else:
      y[ii] = -3*F0/tau_0*t[ii]+3*F0


# Book Problem 3.3
# t = r_[0:8:8000j]       # define the time to look at, easiest to just choose 1 period
# w0 = 2*pi/t[-1]       # define the fundamental frequency (here, I know t(end)=tau)
# tau_0 = 2*pi/w0        # define fundamental period based on w0
# fmax = 1
# y = zeros((len(t),))

# for ii in range(len(t)):
# 	if t[ii] > 1 and t[ii] <= 2:
# 		y[ii] = fmax*(t[ii]-1)
# 	elif t[ii] > 2 and t[ii] <= 3:
# 		y[ii] = -fmax*(t[ii]-2)+fmax
# 	elif t[ii] > 5 and t[ii] <=6:
# 		y[ii] = -fmax*(t[ii]-5)
# 	elif t[ii] > 6 and t[ii] <= 7:
# 		y[ii] = fmax*(t[ii]-6)-fmax


#--------- This is the actual Numerical Fourier Expansion --------
#
# NOTE: SciPy trapz command is trapz(y,t) (opposite order of MATLAB/Octave)

# loop over the number of terms to use in the approximation
def Fourier_Approx(num_terms):
    
    # get the a0 term
    a0 = w0/(2*np.pi)*np.trapz(y,t)  

    # fill arrays with zeros
    a = np.zeros((num_terms,))
    b = np.zeros((num_terms,))
    integral_cos = np.zeros((len(t),num_terms))
    integral_sin = np.zeros((len(t),num_terms))
    sin_term = np.zeros((num_terms,len(t)))
    cos_term = np.zeros((num_terms,len(t)))

    # cycle through the 1 to num_terms Fourier coefficients (a_n and b_n)
    for n in range(num_terms):

        # a_n calculations
        integral_cos[:,n] = y*np.cos((n+1)*w0*t)         # define the integral "interior"
        a[n] = w0/np.pi*np.trapz(integral_cos[:,n],t)   # solve for a_n

        # b_n calculations
        integral_sin[:,n] = y*np.sin((n+1)*w0*t)         # define the integral "interior"
        b[n] = w0/np.pi*np.trapz(integral_sin[:,n],t)   # solve for b_n
    
        sin_term[n,:] = np.sin((n+1)*w0*t)               # calculate the nth sine term
        cos_term[n,:] = np.cos((n+1)*w0*t)               # calculate the nth cosine_term


    # Generate the approximate input based on the Fourier coeff. calculated above
    approx = np.zeros((len(t),)) #First fill with zeros

    for ii in range(len(t)):
         approx[ii] = a0 + np.sum(a*cos_term[:,ii],0) + np.sum(b*sin_term[:,ii],0)

    num_terms
    
    return approx

# Set up basic Figure Properties
#   Many of these setting could also be made default by the .matplotlibrc file
#   axes limits may need to be adjusted depending on input
# fig = plt.figure(figsize=(8,4.5))
fig = plt.figure(figsize=(6,4))
ax = plt.gca()
plt.subplots_adjust(bottom=0.17,left=0.17,top=0.96,right=0.96)
plt.setp(ax.get_ymajorticklabels(),family='CMU Serif',fontsize=18)
plt.setp(ax.get_xmajorticklabels(),family='CMU Serif',fontsize=18)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.grid(True,linestyle=':',color='0.75')
ax.set_axisbelow(True)

plt.xlabel(r'Time (s)',family='CMU Serif',fontsize=22,weight='bold',labelpad=5)
plt.ylabel(r'$f(t)$',family='CMU Serif',fontsize=22,weight='bold',labelpad=10)

exactLine = plt.plot(t,y,'--',linewidth=2,label=r'Exact')

# For the square wave
# plt.ylim(-1.5,2.0)

# For the trapezoid
plt.ylim(0.0,1.5)


f = 'Fourier Approx.'
approxLine, = plt.plot([],[],linewidth=2,label=f)

leg = plt.legend(loc='upper right', fancybox=True)
ltext  = leg.get_texts() 
plt.setp(ltext,family='CMU Serif',fontsize=16)

# place a text box in upper left in axes coords
fourier_text = ax.text(0.04, 0.95, f, transform=ax.transAxes, fontsize=18,
        verticalalignment='top')#, bbox=props)

def init():
    approxLine.set_data([],[])
    fourier_text.set_text('')
    return approxLine,


def animate(i):
    t = np.linspace(0,4,4000)       # define the time to look at, easiest to just choose 1 period
    approx = Fourier_Approx(i+1)
    approxLine.set_data(t,approx)

    f = str(i+1) + '-Terms'

    fourier_text.set_text(f)
#     leg.set_text
#     leg = plt.legend(loc='upper right', fancybox=True)
#     ltext  = leg.get_texts() 
#     print ltext
#     np.shape(ltext)
#     plt.setp(ltext,family='CMU Serif',fontsize=16)
    
    return approxLine, fourier_text,

ani = animation.FuncAnimation(fig, animate, frames=29, init_func=init)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
ani.save('Fourier_Approx.mp4', fps=3)#, extra_args=['-vcodec', 'libx264'])

plt.show()