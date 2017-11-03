# Моделируется сам эксперимент: 1)закачка воды в соседние скважины (сопоставление данных) 2) закачка масла в центр
# сопоставление данных.
import numpy as np
from flow_equations.piezo_in_celindric_2_n_wells import PorePressure_in_Time
from flow_equations.piezo_in_celindric_2_n_wells_PinCenter_with_frac import PorePressure_in_Time as PorePressure_in_Time_2
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
import matplotlib.animation as animation

perm = 200 * 10 ** (-15)  # м2 проницаемость
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
delta_t = 0.01
Pres = 0*10**5
P_center = 65*10**5
q = 0
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3 = k/delta_t
c4 = 1/delta_fi**2
T_exp = 201
T_exp_1 = 20
T_exp_2 = 50
wells_coord = [(int(0.12/delta_r), int(np.pi/2.7/delta_fi)), (int(0.12/delta_r), int((np.pi + np.pi/2.7)/delta_fi))]
P_well = [1500000, 0]
CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений
for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

frac_pressure = [6000000, 5500000, 5000000,4500000, 4000000, 6000000, 5500000, 5000000, 4500000]
frac_angle = 4
frac_angle_2 = 43
frac_coords = [(0, frac_angle), (1, frac_angle), (2, frac_angle), (3, frac_angle),  (0, frac_angle_2), (1, frac_angle_2), (2, frac_angle_2), (3, frac_angle_2), (4, frac_angle_2)]
#frac_coords = []
wells_frac_coords = wells_coord + frac_coords

X = np.zeros((N_r, M_fi))
Y = np.zeros((N_r, M_fi))
for m in range(M_fi):
    for n in range(N_r):
        X[n][m] = (r_well + (n + 1) * delta_r) * np.cos(delta_fi * m)
        Y[n][m] = (r_well + (n + 1) * delta_r) * np.sin(delta_fi * m)

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]


xi = np.linspace(min(X_list), max(X_list), 700)
yi = np.linspace(min(Y_list), max(Y_list), 700)
xig, yig = np.meshgrid(xi, yi)
levels = list(range(0,6500000,10000))

for t in range(T_exp_1):
    print(t)
    Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, mu, perm,
                                        delta_r)
#
#fig = plt.figure()
#P_list = [k for k in Pres_distrib.flat]
#Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')
#surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),
#                            linewidth=0.2, levels=levels)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.show()
#
#perm = 0.05*10**(-15)
#k = mu*fi*(Cf+Cr)/perm
#c3 = k/delta_t
#q = 0.00005
#for t in range(T_exp_2):
#    print(t)
#    Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, mu, perm,
#                                        delta_r)
#    print(min(Pres_distrib.flat), max(Pres_distrib.flat))
#
#fig = plt.figure()
#P_list = [k for k in Pres_distrib.flat]
#Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')
#surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=0, vmax=5000000,
#                            linewidth=0.2, levels=levels)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.show()


def animate(t):
    global Pres_distrib
    global c3,q,perm
    global CP_dict
    print(t)
#    if t < 3:
#        Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, mu, perm,
#                                            delta_r)
#        P_list = [k for k in Pres_distrib.flat]
#        Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')
#
#        surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),
#                            linewidth=0.2, levels=levels)
#
#        return surf,

    if t < 25:
        perm = 0.05 * 10 ** (-15)
        k = mu * fi * (Cf + Cr) / perm
        c3 = k / delta_t
        q = 0.000095
        Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, mu, perm,
                                            delta_r)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))
        P_list = [k for k in Pres_distrib.flat]
        Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')

        surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=1500000,
                            linewidth=0.2, levels=levels)

        return surf,

    if t >= 25:

        for i in range(len(frac_coords)):
            CP_dict[frac_coords[i]] = frac_pressure[i]

        Pres_distrib = PorePressure_in_Time_2(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, delta_r, P_center,
                                              frac_coords, frac_pressure, frac_angle, frac_angle_2, wells_frac_coords)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))
        P_list = [k for k in Pres_distrib.flat]
        Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')

        surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=1500000,
                            linewidth=0.2, levels=levels)

        return surf,

    else:
        return

#def animate(t):
#    global Pres_distrib
#    print(t)
#    Pres_distrib = PorePressure_in_Time_2(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, delta_r, P_center)
#    if t%10 == 0:
#        P_list = [k for k in Pres_distrib.flat]
#        Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')
#
#        surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),
#                            linewidth=0.2, levels=levels)
#
#        return surf,
#    else:
#        return
#
#        #print(min(Pres_distrib.flat), max(Pres_distrib.flat))

fig = plt.figure()
plt.axis('equal')
im_ani = animation.FuncAnimation(fig, animate, 35, interval=1000, blit=False, repeat=False)
#im_ani.save('animation_2st_5exp_with_frac_attempt3.mp4', writer='ffmpeg')
#

plt.show()

