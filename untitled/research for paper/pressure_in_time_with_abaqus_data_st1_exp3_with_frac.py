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

Pinj = 2*10**6
Pprod = 10**5
delta_frac = 0.0005

def pressure_in_abaqus(coord_list_pressure_data):
    x = []
    y = []
    z = []
    with open(coord_list_pressure_data) as coords_file:
        for line in coords_file:
            line = line.strip().split(';')
            if line != []:
                x.append(float(line[0]))
                y.append(float(line[1]))
                z.append(float(line[4]))
            else:
                continue


    xi = np.linspace(min(x),max(x), 500)
    yi = np.linspace(min(y), max(y), 500)
    xig, yig = np.meshgrid(xi, yi)
    zi = interpolate.griddata((x,y), z, (xig, yig), method='cubic')


    return xig, yig, zi

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
    p_in_time = pressure_dict['p'].transpose()*10**6
    Pini = pressure_dict['Pini'].transpose()*10**6
    xhole = pressure_dict['xhole']/1000
    yhole = pressure_dict['yhole']/1000
    print(np.shape(p_in_time), np.shape(x))
    print(list(zip(x[0],y[0])))
    print(xhole,yhole)
    Pinj_row = np.ones((np.shape(p_in_time)[0],1))*Pinj
    Pprod_row = np.ones((np.shape(p_in_time)[0],1))*Pprod

    x = np.append(x, xhole[2])
    x = np.append(x, xhole[0])
    y = np.append(y, yhole[2])
    y = np.append(y, yhole[0])
    Pini = np.append(Pini,Pinj)
    Pini = np.append(Pini, Pprod)
    p_in_time = np.append(p_in_time, Pinj_row, 1)
    p_in_time = np.append(p_in_time, Pprod_row, 1)

    x = np.delete(x, [13]) # убираю нагнетательную скважину, давления в которой определялись Зенченко
    y = np.delete(y, [13])
    p_in_time = np.delete(p_in_time,13,1)
    Pini = np.delete(Pini, [13])
    x = np.delete(x, [12])  # убираю добыв. скважину
    y = np.delete(y, [12])
    p_in_time = np.delete(p_in_time, 12, 1)
    Pini = np.delete(Pini, [12])
    x = np.delete(x, [7])  # убираю датчик с заниженным давлением
    y = np.delete(y, [7])
    p_in_time = np.delete(p_in_time, 7, 1)
    Pini = np.delete(Pini, [7])
    x = np.delete(x, [5]) # убираю датчик с заниженным давлением
    y = np.delete(y, [5])
    p_in_time = np.delete(p_in_time, 5, 1)
    Pini = np.delete(Pini, [5])


    print(np.shape(p_in_time), np.shape(x), np.shape(Pini))
    x = x.reshape(1,np.shape(x)[0])
    y = y.reshape(1,np.shape(y)[0])
#    p_in_time = p_in_time.reshape(np.shape(p_in_time)[0], np.shape(p_in_time)[1])
    print(np.shape(p_in_time), np.shape(x))
    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    xig_2, yig_2 = np.meshgrid(xi, yi)

    return x, y, p_in_time, xig_2, yig_2, Pini

fig2 = plt.figure()
xig,yig,zi = pressure_in_abaqus('PorePr_from_2stExp/PorePr_2st_3exp_2cases.csv')
levels = list(range(0,Pinj+10**5,50000))
surf = plt.contourf(xig, yig, zi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(zi), vmax=np.nanmax(zi),linewidth=0.2, levels = levels)
fig2.colorbar(surf, shrink=0.5, aspect=5)

x, y, p_in_time, xig_2, yig_2, Pini = pressure_in_exp('../Давление/data3.mat')
P_center = p_in_time[:, 9].copy()
print(np.shape(P_center))

print(x[0], y[0], p_in_time[11000])
print(np.shape(x), np.shape(y), np.shape(Pini), np.shape(xig_2))
zi_2 = interpolate.griddata((x[0], y[0]), Pini, (xig_2, yig_2), method='cubic')
plt.contourf(xig_2, yig_2, zi_2, cmap=cm.jet, antialiased=True, vmin=np.nanmin(zi), vmax=np.nanmax(zi), linewidth=0.2, levels = levels)
plt.scatter(x,y)
plt.xlabel('X, m')
plt.ylabel('Y, m')
plt.show()

x = np.delete(x, [9]) # убираю нагнетательную скважину, давления в которой определялись Зенченко
y = np.delete(y, [9])
p_in_time = np.delete(p_in_time, 9, 1)
x = x.reshape(1,np.shape(x)[0])
y = y.reshape(1,np.shape(y)[0])

fig3 = plt.figure()
xig,yig,zi = pressure_in_abaqus('PorePr_from_2stExp/PorePr_2st_3exp_2cases.csv')
levels = list(range(0,Pinj+10**5,50000))
surf = plt.contourf(xig, yig, zi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(zi), vmax=np.nanmax(zi),linewidth=0.2, levels = levels)
fig3.colorbar(surf, shrink=0.5, aspect=5)
k = 0
p_example = list(p_in_time[1])
def animate(i):
    global  k, x, y, p_in_time, p_example
    if 100*i > 11000:
        k += 1
        if k*delta_frac < 0.058:
            x = np.append(x, k*delta_frac*np.cos(np.pi/4 + np.pi/8))
            y = np.append(y, k * delta_frac * np.sin(np.pi / 4 + np.pi / 8))
            p_example.append(P_center[100 * i]/5)
            print(x, p_example)

        if k * delta_frac < 0.049:
            x = np.append(x, k * delta_frac * np.cos(np.pi + np.pi / 4 + np.pi / 8))
            y = np.append(y, k * delta_frac * np.sin(np.pi + np.pi / 4 + np.pi / 8))
            p_example.append(P_center[100 * i]/5)
            print(x, p_example)

    for j in range(np.shape(p_in_time[100*i])[0]):
        p_example[j] = p_in_time[100*i][j]

    p_example_2 = np.array(p_example)
    print(np.shape(p_example_2))
    print(np.shape(x))
    x = x.reshape(1, len(p_example))
    y = y.reshape(1, len(p_example))
    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    xig_2, yig_2 = np.meshgrid(xi, yi)
    print(np.shape(x), np.shape(y), np.shape(p_example_2), np.shape(xig_2))
    zi_2 = interpolate.griddata((x[0], y[0]), p_example_2, (xig_2, yig_2), method='cubic')
    im = plt.contourf(xig_2, yig_2, zi_2, cmap=cm.jet, antialiased=True, vmin=np.nanmin(zi), vmax=np.nanmax(zi), linewidth=0.2, levels = levels)
    print(x, y)
    plt.scatter(x,y)
    print(i)
    return im,

im_ani = animation.FuncAnimation(fig3, animate, 206, interval=1000, blit=False, repeat=False)
#im_ani.save('animation_cubic_interpol_St3Exp2.mp4', writer='ffmpeg')

plt.show()