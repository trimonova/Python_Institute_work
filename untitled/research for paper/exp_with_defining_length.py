import numpy as np
from flow_equations.piezo_in_celindric_2_n_wells import PorePressure_in_Time
from flow_equations.piezo_in_celindric_2_n_wells_with_PinCenter import PorePressure_in_Time as PorePressure_in_Time_2
from flow_equations.piezo_in_celindric_with_wells_frac import PorePressure_in_Time as PorePressure_in_Time_3
from flow_equations.flow_in_frac import Pressure_in_frac
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
import matplotlib.animation as animation



hx = 0.01
hy = 0.01


nu = 0.2
H = 0.07
E = 3*10**9
G = E/2/(1+nu)
k = 4*(1-nu)*H/3.14/G
mu = 0.2
perm = 2*10**(-15)
alpha = 1/24/mu/k/hx**2
#Pinj = 25*10**5
Sh = 5*10**5
Pres = 1*10**5
K1c = 3000

#-----------------------------------------------------------------------------------

mu_water = 100 * 10 ** (-1)  # Па*с вязкость
mu_oil = 100 * 10 ** (-1)
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_oil*fi*(Cf+Cr)/perm
delta_r = 0.01
delta_fi = np.pi / 90 # угол-шаг в радианах
R = 0.215
r_well = 0
N_r = int((R-r_well)/delta_r)
M_fi = int(2*np.pi/delta_fi)
delta_t = 1
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r**2)

q = 0
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3_water = k_water/delta_t
c3_oil = k_oil/delta_t
c3 = k_oil/delta_t
c4 = 1/delta_fi**2
T_exp = 500
T_exp_1 = 1
T_exp_2 = 10
wells_coord = [(int(0.12/delta_r), int(np.pi/2.7/delta_fi)), (int(0.12/delta_r), int((np.pi + np.pi/2.7)/delta_fi))]
P_well = [0, 0]
CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений
for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]
#Pinj = 2000000
#Pprod = 0

#wells_coord = []



def MatFiles_read(FileName):

    f = h5py.File(FileName, 'r')
    pressure_dict = {}
    keys_list = list(f.keys())

    for elem in keys_list:
        k_value = np.array(f.get(elem))
        pressure_dict[elem] = k_value



    print(pressure_dict)

    p_in_time = pressure_dict['p'].transpose() * 10 ** 6



    p_in_time = p_in_time.transpose()


    return p_in_time

frac_angle = int((np.pi/2.7-2*delta_fi)/delta_fi)
frac_angle_2 = int((np.pi+np.pi/2.7 - 7*delta_fi)/delta_fi)
frac_coord_1 = [0,1]
frac_coord_2 = [0,1]
frac_coords = [(int(i /100/ delta_r), frac_angle) for i in frac_coord_1] + [(int(i / 100 / delta_r), frac_angle_2)
                                                                                  for i in frac_coord_2]
l_fr0 = 0.04
l_count = 0
Lfr_mas = [0]
N_fr = int(l_fr0/delta_r)
N_well = 2

FileName = '..\p_hole_3.mat'
P_center = MatFiles_read(FileName)[0]

print(P_center)
fig = plt.figure()
surf = plt.plot(P_center)
plt.show()

q_loss = np.zeros((N_fr - 1, 1))

w = np.zeros((N_fr - 1, 1))
wells_frac_coords = wells_coord + frac_coords

K1_mas = []
for t in range(T_exp):
    w0 = k * (P_center[30000+t*100] - Sh)


    print(t)
    frac_pressure, w = Pressure_in_frac(N_fr, delta_t, N_well, alpha, w0, q_loss, w, k, Sh)

    frac_pressure = np.insert(frac_pressure, N_well, frac_pressure[N_well])
    frac_pressure_new = frac_pressure.copy()
    for i in range(N_well+1):
        frac_pressure_new[i] = frac_pressure[N_well-i]

    #print(frac_pressure)

    coef = -perm / mu * (P_center[10000 + t * 100] / 2 - Pres) / hy / 2
    q_loss = np.ones((N_fr - 1, 1)) * coef
    #q_loss = np.zeros((N_fr - 1, 1))
    for i in range(len(frac_coords)):
        CP_dict[frac_coords[i]] = frac_pressure_new[i]

    #print(frac_coords, frac_pressure)
    Pres_distrib = PorePressure_in_Time_3(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, wells_frac_coords, delta_r, frac_coords, Pres, mu_oil, perm, frac_angle, frac_angle_2, frac_pressure_new, frac_coord_1, frac_coord_2)
    K1 = frac_pressure[int(len(frac_pressure) / 4)] * delta_r / (np.pi * len(frac_pressure) / 2) ** 0.5
    K1_mas.append(K1)
    print(K1)
    if K1 > K1c:
        l_count += 1
        if l_count == 30:
            frac_coord_1.append(frac_coord_1[-1] + 1)
            frac_coord_2.append(frac_coord_2[-1] + 1)
            frac_coords = [(int(i / 100 / delta_r), frac_angle) for i in frac_coord_1] + [
                (int(i / 100 / delta_r), frac_angle_2)
                for i in frac_coord_2]

            w = np.insert(w, 0, 0)
            w = np.insert(w, -1, 0)
            N_fr += 2
            N_well += 1
            wells_frac_coords = wells_coord + frac_coords
            q_loss = np.ones((N_fr - 1, 1)) * coef
            l_count = 0
            Lfr_mas.append(Lfr_mas[-1]+1)
        else:
            Lfr_mas.append(Lfr_mas[-1])
    else:
        Lfr_mas.append(Lfr_mas[-1])




fig = plt.figure()
surf = plt.plot(frac_pressure)
plt.xlabel('Fracture length, cm')
plt.ylabel('Pressure in fracture, Pa')
plt.show()

fig = plt.figure()
surf = plt.plot(Lfr_mas)
plt.axis((0,T_exp, 0, 20))
plt.xlabel('Time, s')
plt.ylabel('Increase in fracture length, *0.5 cm')
plt.show()


X = np.zeros((N_r,M_fi))
Y = np.zeros((N_r, M_fi))
for m in range(M_fi):
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

levels = list(range(0,4000000,10000))
fig = plt.figure()
surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)
plt.xlabel('X, m')
plt.ylabel('Y, m')

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

