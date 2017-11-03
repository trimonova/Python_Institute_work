# Здесь рассматривается прямолинейное поршневое вытеснение маслом воды (Басниев Кочина, стр. 204). Граничные условия: постоянное давление в центральной скважине (левая граница),
# атмосферное давление на правой границе. Значения граничных условий берутся из данных эксперимента. data3-7 - производилась закачка масла преимущественно при постоянном давлении.
# Далее строятся теоритические прямые падения давления от расстояния для двух последовательных временных шагов. При одинаковых граничных условиях. Смотрим, как давление изменилось через некот. время.


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
time = 1000 # s

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

A_2 = delta_p/(mu_oil/mu_water*x_replace - x_replace + L)
A_1 = mu_oil/mu_water*A_2
B_1 = p_center
B_2 = p_bound - A_2*L
x_1 = np.arange(0, x_replace, 0.0001)
x_2 = np.arange(x_replace, L, 0.0001)
p_oil = A_1*x_1 + B_1
p_water = A_2*x_2 + B_2

fig = plt.figure()
plt.plot(x_1, p_oil, 'r')
plt.plot(x_2, p_water, 'b')

const = -(mu_oil-mu_water)*r_well**2/2 - L*mu_water*r_well
D = (mu_water*L/(mu_oil-mu_water))**2 - 2*(perm*delta_p*(time+1000)/fi + const)/(mu_oil-mu_water)
print(D)
a_1 = -L*mu_water/(mu_oil-mu_water) + D**0.5
a_2 = -L*mu_water/(mu_oil-mu_water) - D**0.5
x_replace = max(a_1, a_2)

A_2 = delta_p/(mu_oil/mu_water*x_replace - x_replace + L)
A_1 = mu_oil/mu_water*A_2
B_1 = p_center
B_2 = p_bound - A_2*L
x_1 = np.arange(0, x_replace, 0.0001)
x_2 = np.arange(x_replace, L, 0.0001)

p_oil_2 = A_1*x_1 + B_1
p_water_2 = A_2*x_2 + B_2


plt.plot(x_1, p_oil_2, 'g')
plt.plot(x_2, p_water_2, 'y')

plt.show()
