# 2-D fourier transform of Gaussian function

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def fn(x, y):	#defining the function
	return np.exp(-x*x - y*y)

n_pts = 64		#number of points

x_min = -20		#range of x 
x_max = 20

y_min = -20		#range of y
y_max = 20
dx = (x_max-x_min)/(n_pts - 1)	#spacing
dy = (y_max-y_min)/(n_pts - 1)

x = np.zeros(n_pts)		#declaring the arrays
y = np.zeros(n_pts)		
fxy = np.zeros((n_pts, n_pts))
m = 0

kx = []
ky = []

fk_plot = np.zeros(n_pts*n_pts) 
fk_a = np.zeros(n_pts*n_pts) 

for i in range (0, n_pts):
	x[i] = x_min + i*dx
	for j in range (0, n_pts):
		y[j] = y_min + j*dx
		fxy[i][j] = fn(x[i], y[j])		#setting the values of fn in array		


fk = np.fft.fft2(fxy, norm='ortho')	#numpy fourier transform


k1 = 2*np.pi*(np.fft.fftfreq(n_pts, d = dx))	#k values (x)
k2 = 2*np.pi*(np.fft.fftfreq(n_pts, d = dy))	#k values (y)

for i in range (0, n_pts):
	for j in range (0, n_pts):
		factorx = np.exp(-1j*(k1[i]*x_min))
		factory = np.exp(-1j*(k2[j]*y_min))
		fk_plot[m] = dx*dy*n_pts/(2.0*np.pi)*factorx*factory*fk[i][j]	#fourier transform according to our convention 
		fk_a[m] = 0.5*np.exp(-0.25*(k1[i]*k1[i]+k2[j]*k2[j]))			#analytic fourier transform which is also a Gaussian
		kx.append(k1[i])			#kx axis			
		ky.append(k2[j])			#ky axis
		m = m+1
		

fig = plt.figure()																#3d plot
ax = plt.axes(projection='3d')

ax.scatter(kx, ky, fk_plot, s = 7, cmap='green', label = "Numerical solution")	#plotting numerical fn
ax.scatter(kx, ky, fk_a, s = 5, cmap='green', label = "Analytic solution")		#plotting analytic fn
ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("F(k)")

plt.legend()
plt.savefig("32.png")