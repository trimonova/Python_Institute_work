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
q = 0
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3 = k/delta_t
c4 = 1/delta_fi**2
T_exp = 1500
T_exp_1 = 20
T_exp_2 = 3000
wells_coord = [(int(0.171/delta_r), int(np.pi/4/delta_fi)), (int(0.171/delta_r), int(5*np.pi/4/delta_fi))]
P_well = [2000000, 0]
CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений
for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]



def MatFiles_read(FileName, delta_r, delta_fi):

    f = h5py.File(FileName, 'r')
    pressure_dict = {}
    keys_list = list(f.keys())

    for elem in keys_list:
        k_value = np.array(f.get(elem))
        pressure_dict[elem] = k_value

    x = pressure_dict['xp']/1000
    y = pressure_dict['yp']/1000
    pt = pressure_dict['p']*10**6
    x_list = []
    y_list = []
    for elem in x[0]:
        x_list.append(elem)
    for elem in y[0]:
        y_list.append(elem)
    set_coord = list(zip(x_list, y_list))
    set_with_P = {}
    i = 0
    for couple in set_coord:
        set_with_P[couple] = pt[i]
        i += 1

    set_with_P_mesh = {}
    for key in set_with_P:
        if key != (0, 0):
            r_coord = (key[0]**2+key[1]**2)**0.5

            angle_coord = np.arctan(key[1]/key[0])
            #print(key[0], key[1], r_coord, angle_coord)
            if key[0]<0 and key[1]>=0:
                angle_coord = np.pi + angle_coord
            if key[0]>=0 and key[1]<0:
                angle_coord = 2*np.pi + angle_coord
            if key[0]<0 and key[1]<0:
                angle_coord =np.pi + angle_coord

            r_coord_mesh = round(r_coord/delta_r)
            print(key[0], key[1], r_coord, angle_coord)

            angle_coord_mesh = round(angle_coord/delta_fi)
            set_with_P_mesh[(r_coord_mesh, angle_coord_mesh)] = set_with_P[key]
        else:
            set_with_P_mesh[(0,0)] = set_with_P[(0,0)]

    return set_with_P_mesh



FileName = '..\Давление\data3.mat'
set_with_P_mesh = MatFiles_read(FileName, delta_r, delta_fi)
set_with_P_mesh_theory = {}
for coord in set_with_P_mesh:
    set_with_P_mesh_theory[coord] = []

for t in range(T_exp_1):
    print(t)
    Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, mu, perm,
                                        delta_r, Pres)
#---------------------------------------------------------------------------------------------------------------------

perm = 0.03 * 10 ** (-15)  # м2 проницаемость
k = mu*fi*(Cf+Cr)/perm
c3 = k/delta_t
q = 0.00003

for t in range(T_exp):
    print(t)
    Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, mu, perm, delta_r, Pres)

    for coord in set_with_P_mesh_theory:
        set_with_P_mesh_theory[coord].append(Pres_distrib[coord[0]][coord[1]])

#---------------------------------------------------------------------------------------------------------------------

#frac_pressure = [7500000, 7000000, 6500000,6000000, 5500000, 5000000, 7500000, 7000000, 6500000, 6000000, 5500000]
frac_angle = round((np.pi/4 + np.pi/8)/delta_fi)
frac_angle_2 = round((np.pi + np.pi*3/8)/delta_fi)
frac_coords = [(0, frac_angle), (1, frac_angle), (2, frac_angle), (3, frac_angle),  (4, frac_angle), (5, frac_angle), (0, frac_angle_2), (1, frac_angle_2), (2, frac_angle_2), (3, frac_angle_2), (4, frac_angle_2)]
#frac_coords = []
wells_frac_coords = wells_coord + frac_coords

P_center = 80*10**5

for t in range(T_exp_2):
    print(t)
    P_center = set_with_P_mesh[(0,0)][9500+t]
    frac_pressure = [P_center for i in range(len(frac_coords))]
    for i in range(len(frac_coords)):
        CP_dict[frac_coords[i]] = frac_pressure[i]
    Pres_distrib = PorePressure_in_Time_2(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, delta_r,
                                          P_center,
                                          frac_coords, frac_pressure, frac_angle, frac_angle_2, wells_frac_coords)

    for coord in set_with_P_mesh_theory:
        set_with_P_mesh_theory[coord].append(Pres_distrib[coord[0]][coord[1]])

    #print(min(Pres_distrib.flat), max(Pres_distrib.flat))
time = [i/100 for i in range(4500)]
time1 = [i/100-10 for i in range(4500)]
for coord in set_with_P_mesh:
    fig = plt.figure()

    plt.plot(time, set_with_P_mesh[coord][9500:14000] - set_with_P_mesh[coord][9499])
    plt.plot(time, set_with_P_mesh_theory[coord] - set_with_P_mesh_theory[coord][0], linewidth = 2)
    plt.title((coord[0]*delta_r, coord[1]*delta_fi))
    plt.xlabel('Time, s')
    plt.ylabel('Pressure increase, Pa')
    plt.axis((0, 30, -500,300000))
    plt.show()


