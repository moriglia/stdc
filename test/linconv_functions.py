from stdc.linconv import mult_hermitian, ifft_hermitian, overlap_save_iteration, overlap_save
import numpy as np
from matplotlib import pyplot as plt


N = 64

X = np.empty(shape=(N,), dtype=np.double)
X[0]     = 11
X[2::2] = 3
X[1]     = 5
X[3::2] = 4


H = np.empty_like(X)
H[0]     = 7
H[2::2] = X[2::2]
H[1]     = 13
H[3::2]  = -X[3::2]

print("+------------------------+")
print("| Testing mult_hermitian |")
print("+------------------------+\n")
Y = np.array(mult_hermitian(X, H))

# expected result is
Y_exp = np.empty_like(X)
Y_exp[2::2] = X[2::2]**2 + X[3::2]**2
Y_exp[3::2] = 0
Y_exp[0] = X[0]*H[0]
Y_exp[1] = X[1]*H[1]

print(X[:10:2])
print(X[1:11:2])
print(H[:10:2])
print(H[1:11:2])
print(Y[:10:2])
print(Y[1:11:2])
print(Y_exp[:10:2])
print(Y_exp[1:11:2])
print(f"Result differs from expected result in {sum(Y != Y_exp)} points\n\n\n")

print("+------------------------+")
print("| Testing ifft_hermitian |")
print("+------------------------+\n")

# Preparing the discrete Fourier transform of a decaying exponentian
N_half  = N>>1
theta   = 2*np.pi*np.arange(0,N_half)/N
X[::2]  = (1-0.5*np.cos(theta)) / (5/4 - np.cos(theta))
X[3::2] = (-0.5*np.sin(theta[1:]))  / (5/4 - np.cos(theta[1:]))
X[1]    = 2/3
X *= (1-1/2**N)

x_exp = 0.5**np.arange(0,N)

x = np.array(ifft_hermitian(X))

print(x[0:10])
print(x_exp[0:10])
print(f"Absolute difference of output and expected result is {sum(abs(x-x_exp))}\n\n\n")


print("+--------------------------------+")
print("| Testing overlap_save_iteration |")
print("+--------------------------------+\n")


M = 8

# Fourier discrete transform of the window of length M
theta = np.arange(1, N>>1)*np.pi/N
H[0] = 1
H[1] = np.double(M & 0b1)/M
H[2::2] = np.sin(theta*M)/np.sin(theta)*np.cos(theta*(M-1))
H[3::2] = -np.sin(theta*M)/np.sin(theta)*np.sin(theta*(M-1))

# Input alternating data
x = 1. + (-1)**np.arange(N)

# output should be stuck to the average, i.e. 1
y = np.array(overlap_save_iteration(x, H, M))

print(f"Length of output is {len(y)}, expected is {N-M+1}")

print(x[:10])
print(y[:10])

print(f"Checksum of output: {sum(y)}, which should be {N-M+1}\n\n\n")






print("+----------------------+")
print("| Testing overlap_save |")
print("+----------------------+\n")

Nx = 102_357
x = 1. + (-1)**np.arange(Nx)
h = np.ones(M, dtype=np.double)/M
N = 64

y = np.array(overlap_save(x, h, N))

print(f"Length of convolution is {len(y)}, expected {Nx + M - 1}")



style = {
    "linestyle" : "none",
    "marker" : "x"
}
f  = plt.figure()
subfigs = f.subplots(2,1)
subfigs[0].plot(np.arange(0, N), y[:N] , **style)
subfigs[1].plot(y[-N:Nx], **style)
plt.show()
