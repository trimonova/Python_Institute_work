from scipy import interpolate
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

f = h5py.File(r'..\from_work\data3-7.mat', 'r')

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
i = 0

fig = plt.figure()
#plt.axis([8000, 12000, 0, 10000000])

for dist in dist_from_0:
    if i:
        plt.plot(pt[i])
        i += 1
    else: i += 1
plt.show()

