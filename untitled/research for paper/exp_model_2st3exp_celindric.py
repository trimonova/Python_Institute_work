# Моделируется сам эксперимент: 1)закачка воды в соседние скважины (сопоставление данных) 2) закачка масла в центр
# сопоставление данных.
# пыталась подставлять давления из экспериментов
import numpy as np
from flow_equations.piezo_in_celindric_2_n_wells import PorePressure_in_Time
from flow_equations.piezo_in_celindric_2_n_wells_with_PinCenter import PorePressure_in_Time as PorePressure_in_Time_2
from flow_equations.piezo_in_celindric_with_wells_frac import PorePressure_in_Time as PorePressure_in_Time_3
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
import matplotlib.animation as animation

perm = 200 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2 * 10 ** (-3)
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_oil*fi*(Cf+Cr)/perm
delta_r = 0.01
delta_fi = np.pi / 90 # угол-шаг в радианах
R = 0.215
r_well = 0.0075
N_r = int((R-r_well)/delta_r)
M_fi = int(2*np.pi/delta_fi)
delta_t = 0.01
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r**2)/100
Pres = 0*10**5
q = 0
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3_water = k_water/delta_t
c3_oil = k_oil/delta_t
c3 = k_oil/delta_t
c4 = 1/delta_fi**2
T_exp = 20
T_exp_1 = 10
T_exp_2 = 10
wells_coord = [(int(0.171/delta_r), int(np.pi/4/delta_fi)), (int(0.171/delta_r), int(5*np.pi/4/delta_fi))]
P_well = [2000000, 0]
CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений
for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]
Pinj = 2000000
Pprod = 0

#wells_coord = []



def MatFiles_read(FileName, delta_r, delta_fi):

    f = h5py.File(FileName, 'r')
    pressure_dict = {}
    keys_list = list(f.keys())

    for elem in keys_list:
        k_value = np.array(f.get(elem))
        pressure_dict[elem] = k_value


    pt = pressure_dict['p']*10**6

    x = pressure_dict['xp'] / 1000
    y = pressure_dict['yp'] / 1000
    p_in_time = pressure_dict['p'].transpose() * 10 ** 6
    Pini = pressure_dict['Pini'].transpose() * 10 ** 6
    xhole = pressure_dict['xhole'] / 1000
    yhole = pressure_dict['yhole'] / 1000
    print(np.shape(p_in_time), np.shape(x))
    print(list(zip(x[0], y[0])))
    print(xhole, yhole)
    Pinj_row = np.ones((np.shape(p_in_time)[0], 1)) * Pinj
    Pprod_row = np.ones((np.shape(p_in_time)[0], 1)) * Pprod

    x = np.append(x, xhole[2])
    x = np.append(x, xhole[0])
    y = np.append(y, yhole[2])
    y = np.append(y, yhole[0])
    Pini = np.append(Pini, Pinj)
    Pini = np.append(Pini, Pprod)
    p_in_time = np.append(p_in_time, Pinj_row, 1)
    p_in_time = np.append(p_in_time, Pprod_row, 1)

    x = np.delete(x, [13])  # убираю нагнетательную скважину, давления в которой определялись Зенченко
    y = np.delete(y, [13])
    p_in_time = np.delete(p_in_time, 13, 1)
    Pini = np.delete(Pini, [13])
    x = np.delete(x, [12])  # убираю добыв. скважину
    y = np.delete(y, [12])
    p_in_time = np.delete(p_in_time, 12, 1)
    Pini = np.delete(Pini, [12])
    x = np.delete(x, [7])  # убираю датчик с заниженным давлением
    y = np.delete(y, [7])
    p_in_time = np.delete(p_in_time, 7, 1)
    Pini = np.delete(Pini, [7])
    x = np.delete(x, [5])  # убираю датчик с заниженным давлением
    y = np.delete(y, [5])
    p_in_time = np.delete(p_in_time, 5, 1)
    Pini = np.delete(Pini, [5])
    print(np.shape(p_in_time)[1], np.shape(p_in_time)[0])
    p_in_time = p_in_time.transpose()

    x_list = []
    y_list = []
    for elem in x:
        x_list.append(elem)
    for elem in y:
        y_list.append(elem)
    print(x)
    set_coord = list(zip(x_list, y_list))
    set_with_P = {}
    i = 0
    for couple in set_coord:
        if couple[1] == 0.0 and couple[0] > 0:
            print(couple)
            list_couple = list(couple)
            list_couple[1] = 0.009
            couple = tuple(list_couple)
            print(couple)
        set_with_P[couple] = p_in_time[i]
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

            r_coord_mesh = int(round(r_coord/delta_r))
            #print(key[0], key[1], r_coord, angle_coord)

            angle_coord_mesh = int(round(angle_coord/delta_fi))
            set_with_P_mesh[(r_coord_mesh, angle_coord_mesh)] = set_with_P[key]
        else:
            P_center = set_with_P[(0,0)]

    return set_with_P_mesh, P_center, set_with_P

def graph(Pres_distrib):
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

    xi = np.linspace(min(X_list),max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    print(np.shape(X_list), np.shape(P_list), np.shape(xig))
    Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')
    return xig, yig, Pi

levels = list(range(0,2000000,10000))

for t in range(T_exp):
    print(t)
    #for element in set_with_P_mesh:
    #    CP_dict[element] = set_with_P_mesh[element][t]
    #    print(CP_dict[element])

    Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, mu_water, perm, delta_r, Pres)



xig, yig, Pi = graph(Pres_distrib)
fig = plt.figure()
surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2,
                    levels=levels)

fig.colorbar(surf, shrink=0.5, aspect=5)
#for element in wells_coord:
#    plt.scatter(element[0]*delta_r, element[1]*delta_fi)
plt.xlabel('X, m')
plt.ylabel('Y, m')
plt.show()

wells_coord = []

FileName = '../Давление/data3.mat'
set_with_P_mesh, P_center, set_with_P = MatFiles_read(FileName, delta_r, delta_fi)
for element in set_with_P_mesh:

    wells_coord += [element]
print(wells_coord)
#for t in range(T_exp):
#    print(t)
#    #for element in set_with_P_mesh:
#    #    CP_dict[element] = set_with_P_mesh[element][t]
#    #    print(CP_dict[element])
#    P_center_one = P_center[t]
#
#    Pres_distrib = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, P_center_one, delta_r, delta_fi)
#
#
#
#    xig, yig, Pi = graph(Pres_distrib)
#    fig = plt.figure()
#    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2,
#                        levels=levels)
#
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#    #for element in wells_coord:
#    #    plt.scatter(element[0]*delta_r, element[1]*delta_fi)
#    plt.xlabel('X, m')
#    plt.ylabel('Y, m')
#    plt.show()



#---------------------------------------------------------------------------------------------------------------------
levels = list(range(0,6000000,10000))
mu_water = 2 * 10 ** (-1)  # Па*с вязкость
mu_oil = 2 * 10 ** (-1)

perm = 2 * 10 ** (-15)  # м2 проницаемость
k_water = mu_water*fi*(Cf+Cr)/perm
c3_water = k_water/delta_t

k_oil = mu_oil*fi*(Cf+Cr)/perm
c3_oil = k_oil/delta_t
c3 = k_oil/delta_t
#
#P_center = set_with_P_mesh[(0,0)]
#x_f = 0.05
#N_r_oil = int(x_f/delta_r)
for t in range(T_exp_1):
    print(t)

    for element in set_with_P_mesh:
        CP_dict[element] = set_with_P_mesh[element][9500+t]
        print(CP_dict)
    P_center_one = P_center[9500+t]
    Pres_distrib = PorePressure_in_Time_2(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, P_center_one, delta_r, delta_fi, Pres)


    xig, yig, Pi = graph(Pres_distrib)
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=0, vmax=2000000, linewidth=0.2,
                        levels=levels)
    for element in set_with_P:
        plt.scatter(element[0],element[1])
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.xlabel('X, m')
    plt.ylabel('Y, m')
    plt.show()

#---------------------------------------------------------------------------------------------------------------------
levels = list(range(0,12000000,10000))
#frac_pressure = [7500000, 7000000, 6500000,6000000, 5500000, 5000000, 7500000, 7000000, 6500000, 6000000, 5500000]
frac_angle = round((np.pi/4 + np.pi/8)/delta_fi)
frac_angle_2 = round((np.pi + np.pi*3/8)/delta_fi)
#frac_coord_1 = [0]
#frac_coord_2 = [0]
#frac_coords = [(0, frac_angle), (1, frac_angle), (2, frac_angle), (3, frac_angle),  (4, frac_angle), (5, frac_angle), (0, frac_angle_2), (1, frac_angle_2), (2, frac_angle_2), (3, frac_angle_2), (4, frac_angle_2), (5, frac_angle_2)]
#frac_coords = [(int(i/100/delta_r),frac_angle) for i in frac_coord_1] + [(int(i/100/delta_r),frac_angle_2) for i in frac_coord_2]
#xf = 5
#wells_frac_coords = wells_coord + frac_coords
#print(wells_frac_coords)
q = 0

#perm = 1* 10 ** (-15)  # м2 проницаемость
#k_water = mu_water*fi*(Cf+Cr)/perm
#c3_water = k_water/delta_t
#k_oil = mu_oil*fi*(Cf+Cr)/perm
#c3_oil = k_oil/delta_t

frac_1 = []
frac_2 = []
dlina_1 = 0
dlina_2 = 0
for t in range(T_exp_2):
    print(t)
    P_center_one = P_center[9500+T_exp_1+t]
    if t%2 == 0 and t < 5800:
        frac_1 += [dlina_1]
        dlina_1 += 1
    if t % 2 == 0 and t < 4900:
        frac_2 += [dlina_2]
        dlina_2 += 1
    frac_coord_1 = frac_1
    frac_coord_2 = frac_2
    frac_coords = [(int(i /100/ delta_r), frac_angle) for i in frac_coord_1] + [(int(i / 100 / delta_r), frac_angle_2)
                                                                                  for i in frac_coord_2]
    wells_frac_coords = wells_coord + frac_coords

    frac_pressure_1 = [P_center_one/4 for i in range(len(frac_coord_1))]
    frac_pressure_2 = [P_center_one/4 for i in range(len(frac_coord_2))]
    frac_pressure_1[-1] = P_center_one/4
    frac_pressure_2[-1] = P_center_one/4
    frac_pressure = frac_pressure_1 + frac_pressure_2

    print(wells_frac_coords)
    for i in range(len(frac_coords)):
        CP_dict[frac_coords[i]] = frac_pressure[i]
    print(CP_dict)
    print(np.shape(Pres_distrib))
    Pres_distrib = PorePressure_in_Time_3(N_r, M_fi, Pres_distrib, c1, c2, c3, c4, wells_coord, CP_dict, q, wells_frac_coords, delta_r, frac_coords, Pres, mu_oil, perm, frac_angle, frac_angle_2, frac_pressure)
    print(np.shape(Pres_distrib))

    xig, yig, Pi = graph(Pres_distrib)
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),
                        linewidth=0.2,
                        levels=levels)
    for element in set_with_P:
        plt.scatter(element[0],element[1])

    for element in frac_coords:
        if element[1]*delta_fi < 1.57:
            x_frac = abs(element[0]* np.cos(element[1]*delta_fi)*0.01)
            y_frac = abs(element[0]*np.sin(element[1]*delta_fi)*0.01)
        else:
            x_frac = -abs(element[0] * np.cos(element[1]*delta_fi) * 0.01)
            y_frac = -abs(element[0] * np.sin(element[1]*delta_fi) * 0.01)
        print(element[1],x_frac, y_frac)
        plt.scatter(x_frac, y_frac)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

#fig1 = plt.figure()

#plt.plot(Pres_distrib[:, 30])
#plt.show()

#--------------------------------------------------------------------------------------------------------------
    #print(min(Pres_distrib.flat), max(Pres_distrib.flat))
#time = [i/100 for i in range(4500)]
#time1 = [i/100-10 for i in range(4500)]
#for coord in set_with_P_mesh:
#    fig = plt.figure()
#
#    plt.plot(set_with_P_mesh[coord][9500:12000] - set_with_P_mesh[coord][9499])
#    #plt.plot((set_with_P_mesh_theory[coord] - set_with_P_mesh_theory[coord][0]), linewidth = 2)
#    plt.title((coord[0]*delta_r, coord[1]*delta_fi))
#    plt.xlabel('Time, s')
#    plt.ylabel('Pressure increase, Pa')
#    plt.axis((0, 2500, -10000,400000))
#    plt.show()


