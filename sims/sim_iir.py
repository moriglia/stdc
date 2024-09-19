import numpy as np
from stdc.tdomain import filteriir # filter_butter2_impulseinvariance_iir
import math

from matplotlib import pyplot as plt


D = 10.0 # 10 s
Fs = np.double(5e3) # 5 kHz
N = math.ceil(D*Fs)


f0 = np.double(200) # 200 Hz
w0 = np.double(2.0)*np.pi*f0
w0T = w0/Fs


f1 = np.double(10) # 10 Hz
f2 = np.double(100)
f3 = np.double(200)
f4 = np.double(400)


time_axis = np.linspace(0, N, num=N+1, dtype=np.double)/Fs
tau = 0.5
x = (np.sin(2.0*np.pi*f1*time_axis) + \
     np.sin(2.0*np.pi*f2*time_axis) + \
     np.sin(2.0*np.pi*f3*time_axis) + \
     np.sin(2.0*np.pi*f4*time_axis))* \
     (1-np.exp(-time_axis/tau))

# Normalize power
# x = (N+1)*x/sum(x*x)
Ex = sum(np.abs(x)**2)


a = np.array([1,
              -2.*np.cos(w0T/np.sqrt(2.))*np.exp(-w0T/np.sqrt(2.)),
              np.exp(-np.sqrt(2.)*w0T)],
             dtype=np.double)
b = np.array([0, np.sqrt(2.)*w0T*np.sin(w0T/np.sqrt(2.))*np.exp(-w0T/np.sqrt(2.))],
             dtype=np.double)



y = np.array(filteriir(b, a, x))

# Normalize power:
# y = Ex*y/sum(np.abs(y)**2)


plt.plot(time_axis, x)
plt.plot(time_axis, y)
plt.show()
