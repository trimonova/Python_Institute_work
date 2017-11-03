import numpy as np
import matplotlib.pyplot as plt

perm_water = 2 * 10 ** (-15)  # м2 проницаемость
perm_oil = 2 * 10 ** (-15)
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 0.1
P_0 = 100000
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-7)  # сжимаемость скелета
k_oil = fi*mu_oil*(Cf + Cr)/perm_oil
k_water = fi*mu_water*(Cf + Cr)/perm_water
q = 0.00001
h_oil = 0.0001
h_water = 0.0001
R = 0.215
xf = 0.01
t_step = 1
T = 50
N_oil = int(xf/h_oil)
N_water = int((R-xf)/h_water)

A = np.zeros((N_oil+N_water-1, N_oil+N_water-1))
B = np.zeros((N_oil+N_water-1, 1))
P_old = np.ones((N_oil+N_water-1, 1)) * P_0
print(np.shape(P_old))
for time in range(T):
    for n in range(1, N_oil-1):
        A[n][n - 1] = 1
        A[n][n] = -(2 + k_oil*h_oil**2/t_step)
        A[n][n + 1] = 1
        B[n][0] = -k_oil*h_oil**2/t_step * P_old[n]

    A[0][0] = -(2 + k_oil*h_oil**2/t_step) + 1
    A[0][1] = 1
    B[0][0] = -k_oil*h_oil**2/t_step * P_old[0] - q*mu_oil*h_oil/perm_oil
    A[N_oil - 1][N_oil - 1] = 1/mu_oil/h_oil + 1/mu_water/h_water
    A[N_oil - 1][N_oil - 2] = -1/mu_oil/h_oil
    A[N_oil - 1][N_oil] = -1/mu_water/h_water
    B[N_oil-1][0] = 0

    for m in range(N_oil, N_oil+N_water - 2):
        A[m][m + 1] = 1
        A[m][m - 1] = 1
        A[m][m] = -(2 + k_water * h_water ** 2 / t_step)
        B[m][0] = -k_water * h_water ** 2 / t_step * P_old[n]

    A[N_oil+N_water - 2][N_oil+N_water - 2] = -(2 + k_water*h_water**2/t_step) + 1
    A[N_oil+N_water - 2][N_oil+N_water - 3] = 1
    B[N_oil+N_water-2][0] = -k_water*h_water**2/t_step * P_old[N_oil+N_water-2]

    #print(A)
    #print(B)
    P_new = np.linalg.solve(A,B)
    #print(np.shape(P_new))

    w = (P_new[N_oil-2] - P_new[N_oil-1])/h_oil*perm_oil/mu_oil
    print(w)
    xf_2 = w*t_step*time + xf
    print(xf_2, xf)
    P_old = P_new
    print(P_new)

plt.figure()
plt.plot(P_new)
plt.show()