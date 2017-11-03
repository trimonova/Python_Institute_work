# Решение уравнения пьезопроводности в цилиндрических координатах, двумерный случай.
# в данном файле решается система лин. уравнений. Сначала рассчитывается задача вдоль радиуса, затем вдоль угла.
# В центре скважина с постоянным P, так же задаются координаты n скважин с постоянным давлением.
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
delta_fi = np.pi / 60 # угол-шаг в радианах
R = 0.215
r_well = 0.0075
N_r = int((R-r_well)/delta_r)
M_fi = int(2*np.pi/delta_fi)
delta_t = 0.01
Courant_number = delta_t/k/delta_fi**2 + delta_t/k/delta_r**2
Pres = 1*10**5
P_center = 1000000
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3 = k/delta_t
c4 = 1/delta_fi**2
T_exp = 1
wells_coord = [(int(0.15/delta_r), int(np.pi/4/delta_fi)), (int(0.15/delta_r), int(5*np.pi/4/delta_fi))]
P_well = [1000000, 100000]
CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

frac_pressure = [900000, 800000, 700000,600000, 500000, 400000, 900000, 800000, 700000, 600000, 500000, 400000]
frac_angle = int(3*np.pi/4/delta_fi)
frac_angle_2 = int(np.pi/delta_fi) + frac_angle
frac_coords = [(0, frac_angle), (1, frac_angle), (2, frac_angle), (3, frac_angle), (4, frac_angle), (5, frac_angle), (0, frac_angle_2), (1, frac_angle_2), (2, frac_angle_2), (3, frac_angle_2), (4, frac_angle_2), (5, frac_angle_2)]
#frac_coords = []
wells_frac_coords = wells_coord + frac_coords

for i in range(len(frac_coords)):
    CP_dict[frac_coords[i]] = frac_pressure[i]

def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

#for i in range(len(wells_coord)):
#    Pres_distrib[wells_coord[i][0]][wells_coord[i][1]] = P_well[i]

def PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, delta_r, P_center, frac_coords, frac_pressure, frac_angle, frac_angle_2, wells_frac_coords):

    # пластовое давление во всей области на нулевом временном шаге
    P_total = np.ones((N_r, 1)) * Pres
    for m in range(0, M_fi):

        A = np.zeros((N_r, N_r))
        B = np.zeros((N_r, 1))

        for n in range(1, N_r - 1):
            A[n][n-1] = c1 - c2/((n+1)*delta_r)
            A[n][n] = -2*c1 - c3
            A[n][n+1] = c1 + c2/((n+1)*delta_r)

        A[0][0] = -2*c1 - c3
        A[0][1] = c1 + c2/(1*delta_r)
        A[N_r-1][N_r-1] = -2*c1 - c3 + c1 + c2/((N_r)*delta_r)
        A[N_r-1][N_r-2] = c1 - c2/((N_r)*delta_r)

        if m != M_fi-1:
            for n in range(0, N_r):
                if n == 0:
                    B[n][0] = -c4/((n+1)*delta_r)**2 * Pres_distrib[n][m+1] + 2*c4/((n+1)*delta_r)**2 * Pres_distrib[n][m] - \
                              c4/((n+1)*delta_r)**2 * Pres_distrib[n][m-1] - c3*Pres_distrib[n][m] - \
                              P_center*(c1 - c2/(1*delta_r))

                else:
                    B[n][0] = -c4/((n+1)*delta_r)**2 * Pres_distrib[n][m+1] + 2*c4/((n+1)*delta_r)**2 * Pres_distrib[n][m] - \
                              c4/((n+1)*delta_r)**2 * Pres_distrib[n][m-1] - c3*Pres_distrib[n][m]
        else:
            for n in range(0, N_r):
                if n == 0:
                    B[n][0] = -c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][0] + 2 * c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m] - \
                              c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m - 1] - c3 * Pres_distrib[n][m] - \
                              P_center * (c1 - c2 / (1 * delta_r))

                else:
                    B[n][0] = -c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][0] + 2 * c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m] - \
                              c4 / ((n+1)*delta_r) ** 2 * Pres_distrib[n][m - 1] - c3 * Pres_distrib[n][m]

        wells_coord.sort(key=sortByRad, reverse=True)

        for i in range(len(wells_coord)): # wells_coord начинается с больших радиусов, чтобы удалять с конца
            if m == wells_coord[i][1]:

                A[wells_coord[i][0]+1][wells_coord[i][0]] = 0
                A[wells_coord[i][0]-1][wells_coord[i][0]] = 0
                B[wells_coord[i][0]+1][0] = B[wells_coord[i][0]+1][0] - (c1 - c2/(((wells_coord[i][0]+2)*delta_r)))*CP_dict[wells_coord[i]]
                B[wells_coord[i][0] - 1][0] = B[wells_coord[i][0]-1][0] - (c1 + c2/(((wells_coord[i][0])*delta_r)))*CP_dict[wells_coord[i]]

                A = np.delete(A, wells_coord[i][0], axis=0)
                A = np.delete(A, wells_coord[i][0], axis=1)
                B= np.delete(B, wells_coord[i][0], axis=0)

        half_length = int(len(frac_coords) / 2)

        if m == frac_angle or m == frac_angle_2:
            A[half_length][half_length-1] = 0
            B[half_length][0] = B[half_length][0] - (c1 - c2 / (((half_length + 1) * delta_r))) * \
                                                            frac_pressure[-1]                           # исправила на half_length
            for i in range(half_length-1, -1, -1):
                A = np.delete(A, frac_coords[i][0], axis=0)
                A = np.delete(A, frac_coords[i][0], axis=1)
                B = np.delete(B, frac_coords[i][0], axis=0)

        P_new = np.linalg.solve(A, B)

        wells_coord = wells_coord[::-1]

        if m == frac_angle or m == frac_angle_2:
            for i in range(half_length):
                P_new = np.insert(P_new, frac_coords[i][0], frac_pressure[i])

        for i in range(len(wells_coord)):   # wells_coord начинается с меньших радиусов, чтобы вставлять с начала
            if m == wells_coord[i][1]:
                P_new = np.insert(P_new, wells_coord[i][0], CP_dict[wells_coord[i]])
        #P_new = np.insert(P_new, 0, P_center)
        P_new = P_new.reshape(N_r, 1)

        P_total = np.hstack((P_total, P_new))

    P_total = np.delete(P_total, 0, axis=1)

    Pres_distrib = np.array(P_total.copy())

#-----------------------------------------------------------------------------------

    P_total = np.ones((1, M_fi)) * P_center

    for n in range(0, N_r):

        A = np.zeros((M_fi, M_fi))
        B = np.zeros((M_fi, 1))

        for m in range(1, M_fi-1):
            A[m][m - 1] = c4/((n+1)*delta_r)**2
            A[m][m] = -2*c4/((n+1)*delta_r)**2 - c3
            A[m][m + 1] = c4/((n+1)*delta_r)**2

        A[0][0] = A[m][m]
        A[0][1] = A[m][m+1]
        A[0][M_fi-1] = A[m][m-1]
        A[M_fi-1][M_fi-1] = A[m][m]
        A[M_fi-1][M_fi-2] = A[m][m-1]
        A[M_fi - 1][0] = A[m][m+1]
        #print(A)
        for m in range(M_fi):
           if n == 0:
               B[m][0] = -c1 * (Pres_distrib[n + 1][m] - 2 * Pres_distrib[n][m] + P_center) - \
                         c2 / ((n + 1) * delta_r) * (Pres_distrib[n + 1][m] - P_center) - \
                         c3 * Pres_distrib[n][m]

           elif n == N_r-1:
               B[m][0] = -c1 * (Pres_distrib[n][m] - 2 * Pres_distrib[n][m] + Pres_distrib[n - 1][m]) - \
                         c2 / ((n + 1) * delta_r) * (Pres_distrib[n][m] - Pres_distrib[n - 1][m]) - \
                         c3 * Pres_distrib[n][m]
           else:
               B[m][0] = -c1*(Pres_distrib[n+1][m] - 2*Pres_distrib[n][m] + Pres_distrib[n-1][m]) - \
                         c2/((n+1)*delta_r)*(Pres_distrib[n+1][m] - Pres_distrib[n-1][m]) -\
                         c3*Pres_distrib[n][m]


        #B[0][0] = B[0][0] - c4/((n+1)*delta_r)**2*Pres_distrib[n][M_fi-1]
        #B[M_fi-1][0] = B[M_fi-1][0] - c4/((n+1)*delta_r)**2*Pres_distrib[n][0]

        wells_frac_coords.sort(key=sortByAngle, reverse=True)
        for i in range(len(wells_frac_coords)): # wells_coord начинается с больших углов, чтобы удалять с конца
            if n == wells_frac_coords[i][0]:

                A[wells_frac_coords[i][1]+1][wells_frac_coords[i][1]] = 0
                A[wells_frac_coords[i][1]-1][wells_frac_coords[i][1]] = 0
                B[wells_frac_coords[i][1]+1][0] = B[wells_frac_coords[i][1]+1][0] - c4/((n+1)*delta_r)**2*CP_dict[wells_frac_coords[i]]
                B[wells_frac_coords[i][1] - 1][0] = B[wells_frac_coords[i][1]-1][0] - c4/((n+1)*delta_r)**2*CP_dict[wells_frac_coords[i]]

                A = np.delete(A, wells_frac_coords[i][1], axis=0)
                A = np.delete(A, wells_frac_coords[i][1], axis=1)
                B = np.delete(B, wells_frac_coords[i][1], axis=0)

        #print(B)
        P_new = np.linalg.solve(A, B)

        wells_frac_coords = wells_frac_coords[::-1]
        for i in range(len(wells_frac_coords)):  # wells_coord начинается с меньших углов, чтобы вставлять с начала
            if n == wells_frac_coords[i][0]:

                P_new = np.insert(P_new, wells_frac_coords[i][1], CP_dict[wells_frac_coords[i]])

        P_new = P_new.reshape(M_fi, 1)

        P_total = np.vstack((P_total, P_new.T))

    P_total = np.delete(P_total, 0, axis=0)
    #print(P_total)
    Pres_end = np.array(P_total.copy())

    return Pres_end

if __name__ == '__main__':
    for t in range(T_exp):
        print(t)
        Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord,  CP_dict, delta_r, P_center,  frac_coords, frac_pressure, frac_angle, frac_angle_2, wells_frac_coords)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))

    P_all_center = np.ones((1, M_fi))*P_center
    Pres_distrib = np.vstack((P_all_center, Pres_distrib))
    print(np.shape(Pres_distrib))
    X = np.zeros((N_r+1,M_fi))
    Y = np.zeros((N_r+1, M_fi))
    for m in range(M_fi):
        for n in range(N_r+1):
            X[n][m] = (r_well + (n+1)*delta_r)*np.cos(delta_fi*m)
            Y[n][m] = (r_well + (n+1)*delta_r)*np.sin(delta_fi*m)

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
