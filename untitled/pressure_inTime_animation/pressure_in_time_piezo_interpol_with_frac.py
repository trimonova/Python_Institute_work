import h5py
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy import interpolate

from flow_equations.piezo_for_wells_and_frac_new_BC import PorePressure_in_Time

alpha = 0.8 * 10 ** -12
beta = 0.17 * 10 ** -9
hx = 0.005
hy = 0.005
hz = 0.07

t_step = 1
T_exp = 200
Lx = 0.5
Ly = 0.5

N = int(Lx / hx)  # количество ячеек вдоль оси х
M = int(Ly / hy)
print(N, M)

wells_with_Q = {}

# wells_with_Q = {(int((Lx / 2) / hx), int((Ly / 2 - 0.121) / hy)): -0.00000003}

frac_with_P = {}

Pres = 1 * 10 ** 5  # давление в пласте

V = hx * hy * hz
coeff_1 = hx * hz / hy
coeff_2 = hy * hz / hx
Pres_distrib = np.ones((N, M)) * Pres
Pprod = 100000

f = h5py.File('Pt1.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

print(pressure_dict.keys())

pt = pressure_dict['P']*10**6
t = pressure_dict['t']


g = h5py.File('xpyp.mat', 'r')

coord_dict = {}
keys_list = list(g.keys())

for elem in keys_list:
    k_value = np.array(g.get(elem))
    coord_dict[elem] = k_value

print(coord_dict['yp'])

x = coord_dict['xp']/1000
y = coord_dict['yp']/1000

print(np.shape(y))
print(np.shape(pt))
print(t)

x_list = []
y_list = []

for elem in x[0]:
    x_list.append(elem)

for elem in y[0]:
    y_list.append(elem)

set_coord = list(zip(x_list, y_list))
set_coord.append((0.057, 0.127))
print(set_coord)

set_coord_new_CS = []
angle = 75
angle_rad = np.radians(angle)
print(np.sin(angle_rad))
for couple in set_coord:
    x_new = couple[0]*np.cos(angle_rad) + couple[1]*np.sin(angle_rad) + Lx/2
    y_new = couple[0]*(-np.sin(angle_rad)) + couple[1]*np.cos(angle_rad) + Ly/2
    set_coord_new_CS.append((x_new, y_new))

print(set_coord_new_CS)

def pressure_in_exp(Pprod, t_current, set_coord_new_CS):

    set_with_P = {}
    i = 1
    for couple in set_coord_new_CS:
        try:
            set_with_P[couple] = pt[i][t_current]
            i += 1
        except IndexError:
            set_with_P[couple] = Pprod


    set_with_P_mesh = {}
    for key in set_with_P:
        set_with_P_mesh[(int(key[0] / hx), int(key[1] / hy))] = set_with_P[key]


    return set_with_P_mesh, set_with_P


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
levels = list(range(0, 2000000, 50000))

def animate(t_current):
    set_with_P_mesh, set_with_P = pressure_in_exp(Pprod, 199000+t_current*100, set_coord_new_CS)
    print(t_current)

    global Pres_distrib
    print(np.min(Pres_distrib), np.max(Pres_distrib))
    P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, set_with_P_mesh, frac_with_P, Pres, V, coeff_1, coeff_2, Pres_distrib)
    Pres_distrib = P_total


    P_list = [k for k in P_total.flat]

    Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')
    #print(np.min(P_list), np.max(P_list))


    im = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=0, vmax=2000000,
                        linewidth=0.2, levels=levels)

    for key in set_with_P:
        plt.scatter(key[0], key[1], c='r')


    return im,


   ## surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    #fig.colorbar(surf, shrink=0.5, aspect=5)

fig = plt.figure()

circle = np.arange(0, 2 * np.pi, 0.01)
r = 0.215
plt.plot(r * np.sin(circle) + Lx / 2, r * np.cos(circle) + Ly / 2, 'r')
plt.axis('equal')

im_ani = animation.FuncAnimation(fig, animate, 310, interval=1000, blit=False, repeat=False)
#im_ani.save('animation_st3_exp9_p1_33m10s-38m20s_step_1s.mp4', writer='ffmpeg')

plt.show()



