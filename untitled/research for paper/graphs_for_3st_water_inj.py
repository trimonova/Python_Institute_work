from scipy import interpolate
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

f = h5py.File('Насыщение.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

print(pressure_dict)
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

dist_from_0 = []
for couple in set_coord:
    dist_from_0.append((couple[0]**2 + couple[1]**2)**0.5)

x = []
y = []
i = 1

fig = plt.figure()
plt.axis([800000, 1200000, 0, 1500000])

for dist in dist_from_0:
    if i:
        plt.plot(pt[i])
        i += 1
    else: i += 1
plt.show()

fig2 = plt.figure()
#plt.axis([0, 80, 0, 20])
i = 1
for dist in dist_from_0:
    if i != 11:
        y_min = list(pt[i]).index(np.min(pt[i][6700:7200]))
        x.append(dist*100)
        y.append(y_min/100)
        i += 1
    else: i += 1

print(x)
print(y)
plt.scatter(y, x)


perm = 1.5 * 10 ** (-15)  # м2 проницаемость
mu = 2 * 10 ** (-3)  # Па*с вязкость
B0 = 1.2  # Объемный коэф-т при давление p0
fi = 0.2  # пористость
Cf = 1.6*10 ** (-8)  # сжимаемость флюида
piezo_coeff = perm/mu/B0/fi/Cf
#r = np.arange(0.01, 0.21, 0.01)
#plt.plot(r*100, 1/4/piezo_coeff*r**2)
t_r = np.arange(1, 80, 0.01)
C = 0.01
plt.plot(t_r, (-4*piezo_coeff*t_r*np.log(C*t_r))**0.5*100)

plt.ylabel("Distance to center, cm")
plt.xlabel("Start of pressure increase, s")
plt.show()


fig3 = plt.figure()
#plt.axis([5, 20, 50, 100])
i = 0
for dist in dist_from_0:

    y_min = list(pt[i]).index(np.min(pt[i][7000:11000])) - list(pt[i]).index(np.max(pt[i][7000:11000]))
    x.append(dist*100)
    y.append(-y_min/100)
    i += 1

print(x)
print(y)
plt.scatter(x, y)
r = np.arange(0.01, 0.21, 0.01)
plt.plot(r*100, r*0+69)
plt.ylabel("Time to reach maximum, s")
plt.xlabel("Distance to center, cm")
plt.show()