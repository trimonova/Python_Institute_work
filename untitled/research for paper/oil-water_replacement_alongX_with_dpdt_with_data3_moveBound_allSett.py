import numpy as np
import matplotlib.pyplot as plt
import h5py
import cmath

f = h5py.File(r'..\Давление\data3.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

x = pressure_dict['xp']/1000
y = pressure_dict['yp']/1000
pt = pressure_dict['p']*10**6
xhole = pressure_dict['xhole']/1000
yhole = pressure_dict['yhole']/1000
t = pressure_dict['t']

x_list = []
y_list = []

for elem in x[0]:
    x_list.append(elem)

for elem in y[0]:
    y_list.append(elem)

dist = []

for k in range(len(x_list)):
    dist.append((x_list[k]**2 + y_list[k]**2)**0.5)

set_coord = list(zip(x_list, y_list))
print(set_coord)
#print(dist)


set_with_P = {}
i = 0
for double in set_coord:
    set_with_P[double] = pt[i] # каждой точке с датчиком записываем все данные по давлению за все время эксперимента
    i += 1

perm_water = 2 * 10 ** (-15)  # м2 проницаемость
perm_oil = 4 * 10 ** (-16)
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 0.1

P_0_all = {}
for element in set_coord:
    P_0_all[element] = set_with_P[element][95*100]  # каждой точке с датчиком находим его начальное поровое давление (давление перед закачкой)


frac_end_coords = [(0.058*np.sin(np.pi/16), 0.058*np.cos(np.pi/16)), (-0.049*np.sin(np.pi/16), -0.049*np.cos(np.pi/16)) ] # координаты концов трещины

dist_from_frac_end = [] # расстояние от каждого датчика до ближайшего конца

for element in set_coord:
    dist_1 = ((element[0]-frac_end_coords[0][0])**2 + (element[1]-frac_end_coords[0][1])**2)**0.5
    dist_2 = ((element[0]-frac_end_coords[1][0])**2 + (element[1]-frac_end_coords[1][1])**2)**0.5
    dist_from_frac_end.append(min(dist_1, dist_2))
#print(dist_from_frac_end)
new_dist_from_source = [] # для каждой точки с датчиком выбираем минимальное расстояние до трещины, чтобы определить максимальное влияние закачки на датчик
i = 0
for element in dist:
    new_dist_from_source.append(min(element, dist_from_frac_end[i]))
    print(element, dist_from_frac_end[i], new_dist_from_source[i])
    i += 1

#print(new_dist_from_source)

P_0 = 100000
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-7)  # сжимаемость скелета
k_oil = fi*mu_oil*(Cf + Cr)/perm_oil
k_water = fi*mu_water*(Cf + Cr)/perm_water
q = 0.00004
h_oil = 0.0001
h_water = 0.0001
R = 0.215

t_step = 1
T = 19

xf = 0.0005
N_oil = int(xf/h_oil)
N_water = int((R-xf)/h_water)

A = np.zeros((N_oil+N_water-1, N_oil+N_water-1))
B = np.zeros((N_oil+N_water-1, 1))
P_old = np.ones((N_oil+N_water-1, 1)) * P_0
print(np.shape(P_old))
xf_2 = xf
plt.figure()
for time in range(T):

    N_oil = round(xf / h_oil)
    N_water = round((R - xf) / h_water)
    print(xf, h_oil, xf/h_oil)
    print(N_oil, N_water)
    A = np.zeros((N_oil + N_water - 1, N_oil + N_water - 1))
    B = np.zeros((N_oil + N_water - 1, 1))

    print(np.shape(P_old))

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
        B[m][0] = -k_water * h_water ** 2 / t_step * P_old[m]

    A[N_oil+N_water - 2][N_oil+N_water - 2] = -(2 + k_water*h_water**2/t_step) + 1
    A[N_oil+N_water - 2][N_oil+N_water - 3] = 1
    B[N_oil+N_water-2][0] = -k_water*h_water**2/t_step * P_old[N_oil+N_water-2]

    P_new = np.linalg.solve(A,B)

#    if time == 4:
#        plt.plot(P_new)
#        for i in range(len(dist)):
#
#            plt.scatter(dist[i] * 10000, set_with_P[(x[0][i], y[0][i])][9600 + time * t_step * 100] - P_0_all[(x[0][i], y[0][i])])
#
#    if time == 14:
#        plt.plot(P_new, 'r')
#        for i in range(len(dist)):
#            plt.scatter(dist[i] * 10000, set_with_P[(x[0][i], y[0][i])][9600 + time * t_step * 100] - P_0_all[(x[0][i], y[0][i])])
#        plt.show()

    if time == 18:
        plt.plot(P_new, 'r')
        for i in range(len(dist)):
            plt.scatter(new_dist_from_source[i] * 10000, set_with_P[(x[0][i], y[0][i])][9600 + time * t_step * 100] - P_0_all[(x[0][i], y[0][i])], c='r')
        plt.show()

    w = (P_new[N_oil-2] - P_new[N_oil-1])/h_oil*perm_oil/mu_oil
    print(w)
    xf_2 = w*t_step*time + xf_2
    while xf_2-xf >= h_oil:
        xf = xf + h_oil
    print(xf_2, xf)
    P_old = P_new
    print(P_new)

