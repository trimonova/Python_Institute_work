# Здесь рассматривается прямолинейное поршневое вытеснение маслом воды (Басниев Кочина, стр. 204). Граничные условия: постоянное давление в центральной скважине (левая граница),
# атмосферное давление на правой границе. Значения граничных условий берутся из данных эксперимента. data3-7 - производилась закачка масла преимущественно при постоянном давлении.
# Далее строятся теоритические прямые падения давления от расстояния, Далее на график накладываются значения давлений в точках-датчиках в зависимости от расстояния до центра; отельно строятся графики с данными в точках, лежащих на одной прямой

import numpy as np
import matplotlib.pyplot as plt
import h5py

mu_water = 0.003 # Pa*s
mu_oil = 0.11 # Pa*s
mu_0 = mu_water/mu_oil
L = 0.215 # радиус образца
r_well = 0.0075 # радиус центральной скважины (расстояние от центра, где начиется вода)
perm = 2*10**(-15) # проницаемость
fi = 0.2 # пористость
time = 2000 # s

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

dist = []

for k in range(len(x_list)):
    dist.append((x_list[k]**2 + y_list[k]**2)**0.5)

set_coord = list(zip(x_list, y_list, dist))
print(set_coord)

set_with_P = {}
i = 0
for triple in set_coord:
    set_with_P[triple] = pt[i][time*100]
    i += 1

p_center = set_with_P[0,0,0]
p_bound = 10**5
delta_p = p_bound - p_center

const = -(mu_oil-mu_water)*r_well**2/2 - L*mu_water*r_well
D = (mu_water*L/(mu_oil-mu_water))**2 - 2*(perm*delta_p*time/fi + const)/(mu_oil-mu_water)
print(D)
a_1 = -L*mu_water/(mu_oil-mu_water) + D**0.5
a_2 = -L*mu_water/(mu_oil-mu_water) - D**0.5
x_replace = max(a_1, a_2)
print(x_replace)

print(set_with_P)
A_2 = delta_p/(mu_oil/mu_water*x_replace - x_replace + L)
A_1 = mu_oil/mu_water*A_2
B_1 = p_center
B_2 = p_bound - A_2*L
x_1 = np.arange(0, x_replace, 0.0001)
x_2 = np.arange(x_replace, L, 0.0001)

p_oil = A_1*x_1 + B_1
p_water = A_2*x_2 + B_2

fig = plt.figure()
plt.plot(x_1, p_oil)
plt.plot(x_2, p_water)
for key in set_with_P:
    plt.scatter(key[2], set_with_P[key])
plt.show()

plt.plot(x_1, p_oil)
plt.plot(x_2, p_water)
plt.scatter(0.127, set_with_P[(0,0.127,0.127)])
plt.scatter(0.185, set_with_P[(0,-0.185,0.185)])
plt.scatter(0.07, set_with_P[(0,0.07,0.07)])
plt.scatter(0.07, set_with_P[(0,-0.07,0.07)])
plt.show()

plt.plot(x_1, p_oil)
plt.plot(x_2, p_water)
plt.scatter(0.127, set_with_P[(0.127, 0, 0.127)])
plt.scatter(0.185, set_with_P[(-0.185,0,0.185)])
plt.scatter(0.07, set_with_P[(0.07,0,0.07)])

plt.show()
