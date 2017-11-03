# Решение уравнения пьезопроводности в цилиндрических координатах, двумерный случай.
# в данном файле решается система лин. уравнений только для радиуса - 1ая часть решения двумерного уравнения
# В центре скважина с постоянным q

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm

perm = 2 * 10 ** (-15)  # м2 проницаемость
mu = 2 * 10 ** (-3)  # Па*с вязкость
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k = mu*fi*(Cf+Cr)/perm
delta_r = 0.01
delta_fi = np.pi / 45 # угол-шаг в радианах
R = 0.215
r_well = 0.0075
N_r = int((R-r_well)/delta_r)
M_fi = int(2*np.pi/delta_fi)
delta_t = 0.5
Pres = 1*10**5
q = 0.0001
Pres_distrib = np.ones((N_r, M_fi-1)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3 = k/delta_t
c4 = 1/delta_fi**2
T_exp = 30

print(np.shape(Pres_distrib))
def PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4):
    # пластовое давление во всей области на нулевом временном шаге
    P_total = np.ones((N_r, 1)) * Pres
    for m in range(0, M_fi-1):

        A = np.zeros((N_r, N_r))
        B = np.zeros((N_r, 1))

        for n in range(1, N_r - 1):
            A[n][n-1] = c1 - c2/((n+1)*delta_r)
            A[n][n] = -2*c1 - c3
            A[n][n+1] = c1 + c2/((n+1)*delta_r)

        A[0][0] = -2*c1 - c3 + c1 - c2/(1*delta_r)
        A[0][1] = c1 + c2/(1*delta_r)
        A[N_r-1][N_r-1] = -2*c1 - c3 + c1 + c2/((N_r)*delta_r)
        A[N_r-1][N_r-2] = c1 - c2/((N_r)*delta_r)

        if m != M_fi-2:
            for n in range(0, N_r):
                if n == 0:
                    B[n][0] = -c4/((n+1)*delta_r)**2 * Pres_distrib[n][m+1] + 2*c4/((n+1)*delta_r)**2 * Pres_distrib[n][m] - \
                              c4/((n+1)*delta_r)**2 * Pres_distrib[n][m-1] - c3*Pres_distrib[n][m] - \
                              q*mu*delta_r/perm*(c1 - c2/(1*delta_r))

                else:
                    B[n][0] = -c4/((n+1)*delta_r)**2 * Pres_distrib[n][m+1] + 2*c4/((n+1)*delta_r)**2 * Pres_distrib[n][m] - \
                              c4/((n+1)*delta_r)**2 * Pres_distrib[n][m-1] - c3*Pres_distrib[n][m]
        else:
            for n in range(0, N_r):
                if n == 0:
                    B[n][0] = -c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][0] + 2 * c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m] - \
                              c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m - 1] - c3 * Pres_distrib[n][m] - \
                              q * mu * delta_r / perm * (c1 - c2 / (1 * delta_r))

                else:
                    B[n][0] = -c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][0] + 2 * c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m] - \
                              c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m - 1] - c3 * Pres_distrib[n][m]
        print(B)
        print(A)
        P_new = np.linalg.solve(A, B)

        P_new = P_new.reshape(N_r, 1)

        P_total = np.hstack((P_total, P_new))

    print(np.shape(P_total))
    P_total = np.delete(P_total, 0, axis=1)

    Pres_end = np.array(P_total.copy())
    print(np.shape(Pres_end))
    return Pres_end

if __name__ == '__main__':
    for t in range(T_exp):
        print(t)
        Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))


    X = np.zeros((N_r,M_fi-1))
    Y = np.zeros((N_r, M_fi-1))
    for m in range(M_fi-1):
        for n in range(N_r):
            X[n][m] = n*delta_r*np.cos(delta_fi*m)
            Y[n][m] = n*delta_r*np.sin(delta_fi*m)

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]
    P_list = [k for k in Pres_distrib.flat]


    CP_list = zip(X_list, Y_list, P_list)

    print(min(P_list), max(P_list))

    xi = np.linspace(min(X_list),max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')

    levels = list(range(0,2500000,10000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)

#    t = np.arange(0, 2 * np.pi, 0.01)
#    r = 0.215
#    plt.plot(r * np.sin(t) + Lx/2, r * np.cos(t) + Ly/2)
    #ax = fig.gca(projection='3d')

    #surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
