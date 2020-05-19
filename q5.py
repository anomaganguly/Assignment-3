#Comparing time taken by fft and my method for a range of number of points from 4 to 100
import numpy as np
import time
import matplotlib.pyplot as plt

def fx(x):												#defining fx
	return (np.exp(-x*x))	

nmin = 4												#range of n points
nmax = 100
steps = nmax-nmin + 1



t_mymethod = np.zeros(nmax-nmin + 1)					#arrays for storing time values for different n
t_fft = np.zeros(nmax-nmin + 1)

n_pts = np.linspace(nmin, nmax, num=steps)				#array of n points from 4 to 100

for i in range (0, steps):

	n = nmin + i										#each value of n
	x = np.asarray(np.linspace(-10.0, 10.0, num = n))	#x array
	dx = (x[n-1]-x[0])/(n-1)							#spacing between x points
	kp = 2*np.pi*(np.fft.fftfreq(n, d = dx))			#k array
	
	fk_m = np.zeros(n)									#array to store dft values for my memthod	
	fxarr = np.zeros(n)									#array to store the value of function at x points
	
	for p in range (0, n):
		fxarr[p] = fx(x[p])								#storing function value into an array
										
	t1 = time.time()									#time starts
	for j in range (0, n):	
		factor = 0										#my method starts					
		for l in range (0, n):
			factor = factor + fxarr[l]*np.exp(-1j*kp[j]*x[l])
		fk_m[j] = np.sqrt(1/n)*factor					#dft at each k point
		
	t_mymethod[i] = time.time() - t1					#time taken for my method				
											
	t2 = time.time()									#time starts							
	
	dft = np.fft.fft(fxarr, norm='ortho')				#numpy fft
	
	coeff = np.exp(-1j*kp*x[0])
	
	adft = dx*np.sqrt(n/(2*np.pi))*coeff*dft			#dft with our convention
	t_fft[i] = time.time()-t2							#time taken for fft method



plt.plot(n_pts, t_fft, 'r.', label = "time taken by FFTW")	#plotting
plt.plot(n_pts, t_mymethod, 'b.', label = "time taken by my method")
plt.xlabel("number of points")
plt.ylabel("time (s)")
plt.title("Time taken by different methods")
plt.legend()
plt.show()