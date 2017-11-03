# уравнение пьезопроводности, в одной точке можно задать сток или источник, Q задается вроде бы в м3/с
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

alpha = 0.8*10**-12
beta = 0.17*10**-9
hx = 0.5
hy = 0.5
hz = 10
Q = 0.00001
Pinj = 150*10**5
t_step = 1.5
T_exp = 5000
Lx = 10
Ly = 20
X_well = 5
Y_well = 10
N_well = int(X_well/hx)
M_well = int(Y_well/hy)
N = int(Lx/hx) # количество ячеек вдоль оси х
M = int(Ly/hy)
indic = 0

#wells_with_Q = {int(4/hx,5/hy): -0.0005, int(8/hx, 15/hy): 0.0005}
#wells_with_P = {(int(4/hx),int(5/hy)): 150*10**5, (int(8/hx), int(15/hy)): 80*10**5}
N_well_2 = N_well + 4
M_well_2 = M_well + 4

#print(wells_with_P)

Pres = 100*10**5 # давление в пласте
Pbound = 100*10**5 #  давление на границе


Pres_distrib_0 = np.ones((N+2, M+2))*Pres #  пластовое давление во всей области на нулевом временном шаге
#Pres_distrib_0[N_well_2,M_well_2] = Pinj


V = hx*hy*hz
coeff_1 = hx*hz/hy
coeff_2 = hy*hz/hx


for t in range(1, T_exp):
    P_add_hor = np.ones((N, 1)) * Pres
    P_add_vert = np.ones((1, M + 2)) * Pres
    P_total = np.ones((N, 1)) * Pres
    for m in range(1,M+1):
        A = np.zeros((N,N))
        B = np.zeros((N,1))

        for n in range(1, N-1):
            A[n][n-1] = alpha*coeff_2
            A[n][n] = (-2*coeff_2*alpha - V*beta/t_step)
            A[n][n+1] = alpha*coeff_2

        A[0][0] = -2*coeff_2*alpha - V*beta/t_step
        A[0][1] = alpha*coeff_2
        A[N-1][N-1] = A[0][0]
        A[N-1][N-2] = A[0][1]

        for n in range(1,N+1):
            if n == 1  or n == N:
                B[n-1][0] = -V*beta/t_step*Pres_distrib_0[n][m]- alpha*coeff_1*(Pres_distrib_0[n][m-1] - 2*Pres_distrib_0[n][m] + Pres_distrib_0[n][m+1]) - alpha* coeff_2*Pbound

            elif n == N_well and m == M_well:
                B[n-1][0] = -V*beta/t_step*Pres_distrib_0[n][m]- alpha*coeff_1*(Pres_distrib_0[n][m-1] - 2*Pres_distrib_0[n][m] + Pres_distrib_0[n][m+1]) + Q

            else:
                B[n-1][0] = -V*beta/t_step*Pres_distrib_0[n][m]- alpha*coeff_1*(Pres_distrib_0[n][m-1] - 2*Pres_distrib_0[n][m] + Pres_distrib_0[n][m+1])


        P_new = np.linalg.solve(A,B)
        if indic == 1:
            P_new = np.insert(P_new,N_well_2,Pinj)
            P_new = P_new.reshape(N,1)
            indic = 0

        P_total = np.hstack((P_total,P_new))

    P_total = np.hstack((P_total, P_add_hor))
    P_total = np.vstack((P_add_vert, P_total, P_add_vert))
    Pres_distrib_0 = np.array(P_total.copy())


#----------------------------------------------------------------------------

    P_add_hor = np.ones((N + 2, 1)) * Pres
    P_add_vert = np.ones((1, M)) * Pres
    P_total = np.ones((1, M)) * Pres
    for n in range(1, N + 1):
        A = np.zeros((M, M))
        B = np.zeros((M, 1))

        for m in range(1, M - 1):
            A[m][m - 1] = alpha * coeff_1
            A[m][m] = (-2 * coeff_1 * alpha - V * beta / t_step)
            A[m][m + 1] = alpha * coeff_1

        A[0][0] = -2 * coeff_1 * alpha - V * beta / t_step
        A[0][1] = alpha * coeff_1
        A[M - 1][M - 1] = A[0][0]
        A[M - 1][M - 2] = A[0][1]

        for m in range(1, M + 1):
            if m == 1 or m == M:
                B[m - 1][0] = -V * beta / t_step * Pres_distrib_0[n][m] - alpha * coeff_2 * (Pres_distrib_0[n-1][m] - 2 * Pres_distrib_0[n][m] + Pres_distrib_0[n+1][m]) - alpha * coeff_1 * Pbound

            elif n == N_well and m == M_well:
                B[m - 1][0] = -V * beta / t_step * Pres_distrib_0[n][m] - alpha * coeff_2*(Pres_distrib_0[n-1][ m] - 2 * Pres_distrib_0[n][m] + Pres_distrib_0[n+1][m]) + Q

            else:
                B[m - 1][0] = -V * beta / t_step * Pres_distrib_0[n][m] - alpha * coeff_2 *(Pres_distrib_0[n-1][m] - 2 * Pres_distrib_0[n][m] + Pres_distrib_0[n+1][m])


        P_new = np.linalg.solve(A,B)

        if indic == 1:
            P_new = np.insert(P_new, M_well_2, Pinj)
            P_new = P_new.reshape(M, 1)
            indic = 0

        P_total = np.vstack((P_total, P_new.T))

    P_total = np.vstack((P_total, P_add_vert))
    P_total = np.hstack((P_add_hor, P_total, P_add_hor))

    Pres_distrib_0 = np.array(P_total.copy())

#----------------------------------------------------------------------------

X = np.zeros((N+2,M+2))
Y = np.zeros((N+2, M+2))
for m in range(M+2):
    for n in range(N+2):
        X[n][m] = m*hx
        Y[n][m] = n*hy

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
P_list = [k for k in P_total.flat]

print(X_list)
print(Y_list)
print(P_list)

print(min(P_list), max(P_list))

xi = np.linspace(min(X_list),max(X_list), 700)
yi = np.linspace(min(Y_list), max(Y_list), 700)
xig, yig = np.meshgrid(xi, yi)
Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')


fig = plt.figure()
surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2)
#ax = fig.gca(projection='3d')

#surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()



















