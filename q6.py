# Fourier transform of constant function which chosen to be equal to 1

import numpy as np
import matplotlib.pyplot as plt



n_pts = 128						#number of points

x_min = -50						#range of x
x_max = 50
dx = (x_max-x_min)/(n_pts - 1)	#spacing

x = np.zeros(n_pts)				#declaring the arrays

c = 10
fx = np.dot(c,np.ones(n_pts))	#defining the array of constant fn

k = np.zeros(n_pts)
fk_plot = np.zeros(n_pts) 
fk_a = np.zeros(n_pts) 

fk = np.fft.fft(fx, norm='ortho')#numpy fourier transform


k = 2*np.pi*(np.fft.fftfreq(n_pts, d = dx))		#k values 
factor = np.exp(-1j*k*x_min)
fk_plot = dx*np.sqrt(n_pts/(2*np.pi))*factor*fk	#fourier transform according to our convention

print(fk_plot)
