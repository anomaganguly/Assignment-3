#Convolution of two box functions
import numpy as np
import matplotlib.pyplot as plt

def f(x):			#defining the box fn
	if(abs(x)<=1):
		return 1
	else:
		return 0
		
n_pts = 256			#number of points
x_min = -10			#range of x
x_max = 10
dx = (x_max-x_min)/(n_pts-1)	#spacing

x = np.zeros(n_pts)				#initializing the arrays
f1 = np.zeros(n_pts)
f2 = np.zeros(n_pts)

for i in range (0, n_pts):		#values of function inside an array
	x[i] = x_min + i*dx
	f1[i] = f(x[i])
	f2[i] = f(x[i])	
	
	
fc = dx*np.convolve(f1, f2, 'same')	#convolving the two functions

plt.plot(x, fc, '-o', markersize = 3, label = "convolved function")#plotting convolved fn
plt.plot(x, f1, '-o', markersize = 3, label = "box function")#plotting original fn
plt.xlabel("x")
plt.ylabel("F(x)")
plt.title("Convolution of two box functions")
plt.legend()
plt.show()