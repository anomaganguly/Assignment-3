#Finding DFT and power spectrum of given data 

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

with open('noise.txt') as f:				 #reading the data file
	data = [float(line) for line in f]
	
n_pts = (np.asarray(data).shape)[0]			 #total number of points

dx = 1										 #spacing
x_min = 0

ps = np.zeros(n_pts)						 #declaring array for storing the power spectrum values

x = np.arange(n_pts)						 #x points are 
 
dft = np.fft.fft(data, norm = 'ortho')		 #numpy discrete fourier transform
k = 2*np.pi*(np.fft.fftfreq(n_pts, d = dx))	 #k values 

factor = np.exp(-1j*k*x_min)
adft = dx*np.sqrt(n_pts/(2*np.pi))*factor*dft #dft according to our convention


for i in range (0, n_pts):
	ps[i] = (1/n_pts)*adft[i]*np.conjugate(adft[i]) #power spectrum

#plt.plot(x, data, 'b.', markersize = 3)
#plt.plot(k, adft, 'b.', markersize = 3)
#plt.title("Data set")
plt.title("Histogram of power spectrum of given data")

#plt.plot(k, ps, 'b.', markersize = 3)

plt.hist(ps, bins=10)
plt.show()