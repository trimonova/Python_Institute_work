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
mu = 0.11 # Па*с вязкость
B0 = 1.2  # Объемный коэф-т при давление p0
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 10 ** (-10)  # сжимаемость скелета
B = B0  # на самом деле: B = B0/(1+Cf*(p-p0)) - объемный коэф-т при давлении p
fi0 = fi  # пористость при давлении p0
alpha = perm / mu / B
beta = fi * Cf / B0 + fi0 * Cr / B
eta = alpha / beta


hx = 0.005
hy = 0.005
hz = 0.07

t_step = 0.01
T_exp = 70
Lx = 0.5
Ly = 0.5

N = int(Lx / hx)  # количество ячеек вдоль оси х
M = int(Ly / hy)
print(N, M)

wells_with_Q = {(N/2, M/2): -0.37*10**(-8)}
wells_with_P_mesh = {}

# wells_with_Q = {(int((Lx / 2) / hx), int((Ly / 2 - 0.121) / hy)): -0.00000003}

frac_with_P = {}


V = hx * hy * hz
coeff_1 = hx * hz / hy
coeff_2 = hy * hz / hx

Pprod = 0
Pinj = 10*10**5

Pres = (Pinj-Pprod)/2  # давление в пласте
Pres_distrib = np.ones((N, M)) * Pres

f = h5py.File('../Давление/data1.mat', 'r')

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
Pini = pressure_dict['Pini'].transpose()*10**6

x_list = []
y_list = []

for elem in x[0]:
    x_list.append(elem)

for elem in y[0]:
    y_list.append(elem)

set_coord = list(zip(x_list, y_list))
print(set_coord)

center_location = set_coord.index((0,0))
print(center_location)

X = np.zeros((N, M))
Y = np.zeros((N, M))
for m in range(M):
    for n in range(N):
        X[n][m] = n * hx
        Y[n][m] = m * hy

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
xi = np.linspace(min(X_list), max(X_list), 1000)
yi = np.linspace(min(Y_list), max(Y_list), 1000)
xig, yig = np.meshgrid(xi, yi)
levels = list(range(0, Pinj+10**5, 50000))
P_in_central_point_list = []

def animate(t_current):

    global Pres_distrib
    global P_in_central_point_list

    P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, wells_with_P_mesh, frac_with_P, Pres, V, coeff_1, coeff_2, Pres_distrib)
    Pres_distrib = P_total

    P_in_central_point = P_total[N/2][M/2]

    P_in_central_point_list.append(P_total[N/2][M/2])

    #print('P_in_central_point ', P_in_central_point)
    #print('P_exp ', pt[center_location][t_current + 6800])

    if t_current%100 == 0:
        average_value = sum(P_in_central_point_list)/len(P_in_central_point_list)
        max_value = np.max(pt[center_location][t_current+6800:t_current+6800+100])
        print ('average_Theory ', average_value)
        print('max_exp ', max_value)
        P_in_central_point_list = []

    print(t_current)

    P_list = [k for k in P_total.flat]

    Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')

    im = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=0, vmax=1100000,
                        linewidth=0.2, levels=levels)



    return im,


fig = plt.figure()

circle = np.arange(0, 2 * np.pi, 0.01)
r = 0.215
plt.plot(r * np.sin(circle) + Lx / 2, r * np.cos(circle) + Ly / 2, 'r')
plt.axis('equal')

im_ani = animation.FuncAnimation(fig, animate, 3500, interval=100, blit=False, repeat=False)
#im_ani.save('animation_st2_attempt2.mp4', writer='ffmpeg')

plt.show()