#! /usr/bin/env python 

##########################################################################################
# Script to complete a numerical approximation of a Fourier Expansion
#
#                         ----- THIS IS NOT AN FFT!!! -----
# 
# Created: 3/14/13 
#   - Joshua Vaughan 
#   - joshua.vaughan@louisiana.edu
#   - http://www.ucs.louisiana.edu/~jev9637
#
# Modified:
#   * 3/14/13 - Joshua Vaughan - joshua.vaughan@louisina.edu)
#       - better commenting
#       - added a section to compare the response of a mass-spring-damper to the exact 
#           input to the reponse to the approximated input, it is commented out by default
#
##########################################################################################

from numpy import *             # Grab all of the NumPy functions
from matplotlib.pyplot import * # Grab MATLAB plotting functions
import control                  # import the control system functions


#--------- Input your function to examine here --------
#

# This is the square wave from the Fall 2012 midterm example we did in class
t = r_[0:4:8000j]       # define the time to look at, easiest to just choose 1 period
w0 = 2*pi/t[-1]       # define the fundamental frequency (here, I know t(end)=tau)
tau_0 = 2*pi/w0        # define fundamental period based on w0
#
y = (t>-4) - 2*(t>-2) + 2*(t>0) - 2*(t>2) + 2*(t>4) - 2*(t>6) + 2*(t>8) - (t>10)

# This is the "trapezoid" wave example we worked in class
# t = r_[0:4:4000j]       # define the time to look at, easiest to just choose 1 period
# w0 = 2*pi/t[-1]       # define the fundamental frequency (here, I know t(end)=tau)
# tau_0 = 2*pi/w0        # define fundamental period based on w0
#
# F0 = 1
# y = zeros((len(t),))
# 
# for ii in range(len(t)):
#   if t[ii] <= tau_0/3:
#       y[ii] = 3*F0/tau_0*t[ii]
#   elif t[ii] <= 2*tau_0/3:
#       y[ii] = F0
#   else:
#       y[ii] = -3*F0/tau_0*t[ii]+3*F0


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

# define the number of terms to use in the approximation
num_terms = 3  

# get the a0 term
a0 = w0/(2*pi)*trapz(y,t)  

# fill arrays with zeros
a = zeros((num_terms,))
b = zeros((num_terms,))
integral_cos = zeros((len(t),num_terms))
integral_sin = zeros((len(t),num_terms))
sin_term = zeros((num_terms,len(t)))
cos_term = zeros((num_terms,len(t)))

# cycle through the 1 to num_terms Fourier coefficients (a_n and b_n)
for n in range(num_terms):

    # a_n calculations
    integral_cos[:,n] = y*cos((n+1)*w0*t)         # define the integral "interior"
    a[n] = w0/pi*trapz(integral_cos[:,n],t)   # solve for a_n

    # b_n calculations
    integral_sin[:,n] = y*sin((n+1)*w0*t)         # define the integral "interior"
    b[n] = w0/pi*trapz(integral_sin[:,n],t)   # solve for b_n
    
    sin_term[n,:] = sin((n+1)*w0*t)               # calculate the nth sine term
    cos_term[n,:] = cos((n+1)*w0*t)               # calculate the nth cosine_term


# Generate the approximate input based on the Fourier coeff. calculated above
approx = zeros((len(t),)) #First fill with zeros

for ii in range(len(t)):
     approx[ii] = a0 + sum(a*cos_term[:,ii],0) + sum(b*sin_term[:,ii],0)

#   Many of these setting could also be made default by the .matplotlibrc file
#   axes limits may need to be adjusted depending on input
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

xlabel(r'Time (s)',family='CMU Serif',fontsize=22,weight='bold',labelpad=5)
ylabel(r'f(t)',family='CMU Serif',fontsize=22,weight='bold',labelpad=10)

plot(t,y,'--',linewidth=2,label=r'Exact')

f = str(num_terms) + '-Term Fourier Expansion'
plot(t,approx,linewidth=2,label=f)

leg = legend(loc='upper right', fancybox=True)
ltext  = leg.get_texts() 
setp(ltext,family='CMU Serif',fontsize=16)

# save the figure as a high-res pdf in the current folder
f = str(num_terms) + 'orderFourierApprox.pdf'
savefig(f,dpi=600)




#--------- Uncomment below to compare the response to the exact function with --------
#---------  our n-term fourier approximation --------
#
wn = (2*pi)**2
z = 0.1

# Define the system to use in simulation - in transfer function form here
num = [2.*z*wn,wn**2];
den = [1,2.*z*wn,wn**2];

sys = control.tf(num,den);
       
# run the simulation - first with the exact input
[T_out,yout_exact,xout_exact] = control.forced_response(sys,t,y)

# run the simulation - utilize the built-in initial condition response function
[T_approx,yout_approx,xout_approx] = control.forced_response(sys,t,approx)

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
ylabel('Response',family='CMU Serif',fontsize=22,weight='bold',labelpad=10)

plot(t,yout_exact,'--',linewidth=2,label=r'Exact')

f = str(num_terms) + '-Term Fourier Approx. Resp.'
plot(t,yout_approx,linewidth=2,label=f)

# ylim(-1.25*min(yout_exact),1.5*max(yout_exact))

leg = legend(loc='upper right', fancybox=True)
ltext  = leg.get_texts() 
setp(ltext,family='CMU Serif',fontsize=16)

# save the figure as a high-res pdf in the current folder
f = str(num_terms) + 'orderFourierApprox_Resp.pdf'
savefig(f,dpi=600)

# show the figure
show()


