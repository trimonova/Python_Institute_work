from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import types
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

wells_with_Q = {}
frac_with_P = {}
Pres = 0 * 10 ** 5  # давление в пласте
Pprod = 0
Pinj = 1450000
hx = 0.005
hy = 0.005
hz = 0.07
t_step = 0.05
T_exp = 70
Lx = 0.5
Ly = 0.5

N = int(Lx / hx)  # количество ячеек вдоль оси х
M = int(Ly / hy)

V = hx * hy * hz
coeff_1 = hx * hz / hy
coeff_2 = hy * hz / hx
Pres_distrib = np.ones((N, M)) * Pres
wells_with_P_mesh = {}
wells_with_P_mesh[(int((0.057 + Lx / 2) / hx), int((0.127 + Ly / 2) / hy))] = 1450000
wells_with_P_mesh[(int((-0.057 + Lx / 2) / hx), int((-0.127 + Ly / 2) / hy))] = 0

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
levels = list(range(0, 1600000, 50000))


for times in range(T_exp):
    t_step = 0.5
    P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, wells_with_P_mesh, frac_with_P, Pres, V,
                                   coeff_1, coeff_2, Pres_distrib)
    Pres_distrib = P_total

P_list = [k for k in P_total.flat]
Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')

def pressure_in_exp(pore_pressure_mat):

    f = h5py.File(pore_pressure_mat, 'r')

    pressure_dict = {}
    keys_list = list(f.keys())

    for elem in keys_list:
        k_value = np.array(f.get(elem))
        pressure_dict[elem] = k_value
    print(pressure_dict.keys())
    x = pressure_dict['xp']/1000
    y = pressure_dict['yp']/1000
    p_in_time = pressure_dict['P'].transpose()*10**6

    xhole = pressure_dict['xhole']/1000
    yhole = pressure_dict['yhole']/1000
    print(np.shape(p_in_time), np.shape(x))
    print(list(zip(x[0],y[0])))
    print(xhole,yhole)
    Pinj_row = np.ones((np.shape(p_in_time)[0],1))*Pinj
    Pprod_row = np.ones((np.shape(p_in_time)[0],1))*Pprod

    x = np.append(x, xhole[2]) # добавляю нагнетат. и доб. скважины
    x = np.append(x, xhole[0])
    y = np.append(y, yhole[2])
    y = np.append(y, yhole[0])

    p_in_time = np.append(p_in_time, Pinj_row, 1)
    p_in_time = np.append(p_in_time, Pprod_row, 1)
    p_in_time = np.delete(p_in_time, 0, 1)

    print(np.shape(p_in_time), np.shape(x))
    x = x.reshape(1,np.shape(x)[0])+Lx/2
    y = y.reshape(1,np.shape(y)[0])+Ly/2
#    p_in_time = p_in_time.reshape(np.shape(p_in_time)[0], np.shape(p_in_time)[1])
    print(np.shape(p_in_time), np.shape(x))
    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    xig_2, yig_2 = np.meshgrid(xi, yi)

    return x, y, p_in_time, xig_2, yig_2

fig2 = plt.figure()

levels = list(range(0,Pinj+10**5,50000))
surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels = levels)
fig2.colorbar(surf, shrink=0.5, aspect=5)

x, y, p_in_time, xig_2, yig_2 = pressure_in_exp('Насыщение.mat')


def animate(i):
    print(np.shape(x), np.shape(y),  np.shape(p_in_time))
    zi_2 = interpolate.griddata((x[0], y[0]), p_in_time[100*i], (xig_2, yig_2), method='cubic')
    im = plt.contourf(xig_2, yig_2, zi_2, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2, levels = levels)
    plt.scatter(x,y)
    print(i)
    return im,

im_ani = animation.FuncAnimation(fig2, animate, 206, interval=1000, blit=False, repeat=False)
#im_ani.save('animation_cubic_interpol_St3Exp2.mp4', writer='ffmpeg')

plt.show()