import numpy as np
import matplotlib.pyplot as plt
import sympy


perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 0.11
P_0 = 100000
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_oil = fi*mu_oil*(Cf + Cr)/perm
k_water = fi*mu_water*(Cf + Cr)/perm
q = 0.000001

R = 0.215
xf = sympy.Symbol('xf')
t_step = 0.5
N_oil = 3
N_water = 3
h_oil = xf/N_oil
h_water = (R-xf)/N_water
print(h_oil, h_water)

A = sympy.zeros(N_oil+N_water-1)
B = sympy.zeros(N_oil+N_water-1, 1)
P_old = sympy.ones(N_oil+N_water-1, 1) * P_0
print(type(A))

for n in range(1, N_oil-1):
    A[n, n-1] = 1
    A[n, n] = -(2 + k_oil*h_oil**2/t_step)
    A[n, n + 1] = 1
    B[n, 0] = -k_oil*h_oil**2/t_step * P_old[n]

print(A)
A[0, 0] = -(2 + k_oil*h_oil**2/t_step) + 1
A[0, 1] = 1
B[0, 0] = -k_oil*h_oil**2/t_step * P_old[0] - q*mu_oil*h_oil/perm
A[N_oil - 1, N_oil - 1] = 1/mu_oil/h_oil + 1/mu_water/h_water
A[N_oil - 1, N_oil - 2] = -1/mu_oil/h_oil
A[N_oil - 1, N_oil] = -1/mu_water/h_water
B[N_oil-1, 0] = 0

for m in range(N_oil, N_oil+N_water - 2):
    A[m, m + 1] = 1
    A[m, m - 1] = 1
    A[m, m] = -(2 + k_water * h_water ** 2 / t_step)
    B[m, 0] = -k_water * h_water ** 2 / t_step * P_old[n]

A[N_oil+N_water - 2, N_oil+N_water - 2] = -(2 + k_water*h_water**2/t_step) + 1
A[N_oil+N_water - 2, N_oil+N_water - 3] = 1
B[N_oil+N_water-2, 0] = -k_water*h_water**2/t_step * P_old[N_oil+N_water-2]

print(A)
print(B)
A_inv = A.inv()
print(A_inv)
P_new = A_inv*B
print(P_new)
xf_new = (P_new[N_oil-2] - P_new[N_oil-1])/h_oil*perm/mu_oil*t_step

print(sympy.limit(xf_new, xf, 5*10**(-7)))

#podintegral_equation = 1/dif
#xf_new = sympy.Symbol('xf_new')
#integral_equation = sympy.integrate(1/dif, (xf, 0, xf_new)) - t_step
#
#print(integral_equation)

#P_new = np.linalg.solve(A,B)
