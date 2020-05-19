# Fourier transform of sinc function using NumPy

import numpy as np
import matplotlib.pyplot as plt

def fn(x):						#defining the function
	if(x != 0):
		return np.sin(x)/x
	else:
		return 1


n_pts = 256						#number of points

x_min = -50						#range of x
x_max = 50
dx = (x_max-x_min)/(n_pts - 1)	#spacing

x = np.zeros(n_pts)				#declaring the arrays
x[0] = x_min
fx = np.zeros(n_pts)
k = np.zeros(n_pts)
fk_plot = np.zeros(n_pts) 
fk_a = np.zeros(n_pts) 

for i in range (0, n_pts):
	fx[i] = fn(x_min + i*dx)	#setting the values of fn in array
	x[i] = x_min + i*dx


fk = np.fft.fft(fx, norm='ortho')#numpy fourier transform


k = 2*np.pi*(np.fft.fftfreq(n_pts, d = dx))		#k values 
factor = np.exp(-1j*k*x_min)
fk_plot = dx*np.sqrt(n_pts/(2*np.pi))*factor*fk	#fourier transform according to our convention


for i in range (0, n_pts):						#analytic fourier transform which is a box fn
	if (abs(k[i])<=1):
		fk_a[i] = np.sqrt(2*np.pi)/2
	else:
		fk_a[i] = 0

################################################################################################

file = open("fftw_1.dat", "r")					#reading the file which has fourier transform done using fftw  
lines = file.readlines()
file.close()
kfftw = []
fftw = []

for line in lines[1:]:
	p = line.split()
	kfftw.append(float(p[0]))
	fftw.append(float(p[1]))
	
kfftw = np.asarray(kfftw)
fftw = np.asarray(fftw)

####################################################################################################

file = open("fftgsl.dat", "r")					#reading the file which has fourier transform done using gsl 
lines = file.readlines()
file.close()
kgsl = []
gsl = []

for line in lines[1:]:
	q = line.split()
	kgsl.append(float(q[0]))
	gsl.append(float(q[1]))
	
kgsl = np.asarray(kgsl)
gsl = np.asarray(gsl)

########################################################################################################

plt.plot(k, fk_plot, '-o', linewidth = 0.5, markersize = 3, label = "NumPy solution")	#plotting
plt.plot(k, fk_a, '-o', linewidth = 0.5, markersize = 3, label = "Analytic solution")
plt.plot(kfftw, fftw, '-o', linewidth = 0.5, markersize = 3, label = "FFTW solution")
plt.plot(kgsl, gsl, '-o', linewidth = 0.5, markersize = 3, label = "GSL solution")
plt.xlabel("k")
plt.ylabel("F(k)")
plt.title("Fourier transform of sinc function")
plt.legend()
plt.show()

