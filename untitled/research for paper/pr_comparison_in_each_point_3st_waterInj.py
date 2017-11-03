from scipy import interpolate
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

from flow_equations.piezo_for_wells_and_frac_new_BC import PorePressure_in_Time

perm = 2 * 10 ** (-15)  # м2 проницаемость
mu = 2 * 10 ** (-3)  # Па*с вязкость
B0 = 1.2  # Объемный коэф-т при давление p0
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
B = B0  # на самом деле: B = B0/(1+Cf*(p-p0)) - объемный коэф-т при давлении p
fi0 = fi  # пористость при давлении p0
alpha = perm / mu / B
beta = fi * Cf / B0 + fi0 * Cr / B
eta = alpha / beta


hx = 0.005
hy = 0.005
hz = 0.07


T_exp = 2000
Lx = 0.5
Ly = 0.5

N = int(Lx / hx)  # количество ячеек вдоль оси х
M = int(Ly / hy)
print(N, M)

wells_with_Q = {}

# wells_with_Q = {(int((Lx / 2) / hx), int((Ly / 2 - 0.121) / hy)): -0.00000003}

frac_with_P = {}

Pres = 0 * 10 ** 5  # давление в пласте

V = hx * hy * hz
coeff_1 = hx * hz / hy
coeff_2 = hy * hz / hx
Pres_distrib = np.ones((N, M)) * Pres
Pprod = 0
Pinj = 14.5*10**5

f = h5py.File('Насыщение.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

x = pressure_dict['xp']/1000
y = pressure_dict['yp']/1000
pt = pressure_dict['P']*10**6
xhole = pressure_dict['xhole']/1000
yhole = pressure_dict['yhole']/1000
t = pressure_dict['t']
print(np.shape(y))
print(np.shape(pt))

x_list = []
y_list = []

for elem in x[0]:
    x_list.append(elem)

for elem in y[0]:
    y_list.append(elem)

set_coord = list(zip(x_list, y_list))
print(set_coord)

set_with_P = {}
i = 1
for couple in set_coord:
    set_with_P[couple] = pt[i]
    i += 1

set_with_P_mesh = {}
for key in set_with_P:
    set_with_P_mesh[(int((key[0] + Lx / 2) / hx), int((key[1] + Ly / 2) / hy))] = set_with_P[key]


wells_with_P = {}
wells_with_P[(xhole[0][0], yhole[0][0])] = Pprod
wells_with_P[(xhole[2][0], yhole[2][0])] = Pinj
wells_with_P_mesh = {}
for key in wells_with_P:
    wells_with_P_mesh[(int((key[0]+Lx/2)/hx), int((key[1]+Ly/2)/hy))] = wells_with_P[key]

P_in_setting_points = {}
for key in set_with_P_mesh:
    P_in_setting_points[(key[0], key[1])] = []

print(P_in_setting_points)

for times in range(T_exp):
    t_step = 0.1
    P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, wells_with_P_mesh, frac_with_P, Pres, V,
                                   coeff_1, coeff_2, Pres_distrib)
    Pres_distrib = P_total

    for key in set_with_P_mesh:
        P_in_setting_points[(key[0],key[1])].append(P_total[key[0]][key[1]])


time_1 = np.arange(0, 1330204)/100
time_2 = np.arange(0, 2000)/10

levels = list(range(0, 1600000, 50000))

for key in set_with_P_mesh:
    print(key)
    fig = plt.figure()
    plt.axis([0, 500, 0, 1300000])
    plt.plot(time_1, set_with_P_mesh[key])
    plt.plot(time_2, P_in_setting_points[key])


    plt.show()



