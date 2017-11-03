# Здесь рассматривается прямолинейное поршневое вытеснение маслом воды (Извеков + Басниев Кочина, стр. 204).
#За Нестационарность вычислений отвечает меняющиеся граничное условие на скважине. Рассматривается почти аналитическое решение, кроме части, отвечающей за изменение давления на забое.
# учитывается изменение давление на скважине.
# атмосферное давление на правой границе. Значения граничных условий берутся из данных эксперимента. data3-7 - производилась закачка масла преимущественно при постоянном давлении.
# строятся графики, количество которых равно кол. точек-датчиков. Берутся экспериментальные данные (изменение давления во времени) в каждой точке. Строится график.
# на Этот же график накладывается теоретическая кривая, соответствующая изменению давления в точке во времени.
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

t_step = 1
Texp = int(1500/t_step)

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

set_coord = list(zip(x_list, y_list))
print(set_coord)

set_with_P = {}
i = 0
for double in set_coord:
    set_with_P[double] = pt[i]
    i += 1

p_center = set_with_P[0,0]
p_bound = 10**5

x_replace = 0.0075

delta_p_2 = p_bound - p_center[t_step*100]
integral = delta_p_2/2*t_step
const = -perm/fi*integral - L*mu_water*x_replace - (mu_oil - mu_water)*x_replace**2/2

x_dist_pres = {}
for i in range(len(dist)):
    x_dist_pres[(set_coord[i], dist[i])] = []

x_dist_pres_exp = {}
for i in range(len(dist)):
    x_dist_pres_exp[(set_coord[i], dist[i])] = []

print(x_dist_pres_exp)


for time in range(1,Texp+1):

    delta_p_1 = p_bound - p_center[time*100*t_step]
    delta_p_2 = p_bound - p_center[(time+1)*100*t_step]
    integral += (delta_p_1+delta_p_2)/2*t_step
    D = (L*mu_water/(mu_oil - mu_water))**2 - 2*const/(mu_oil-mu_water) - 2*perm/fi/(mu_oil-mu_water)*integral
    a_1 = (-L*mu_water/(mu_oil - mu_water) + D**0.5)
    a_2 = (-L*mu_water/(mu_oil - mu_water) - D**0.5)

    x_replace = max(a_2, a_1)

    A_2 = delta_p_1 / (mu_oil / mu_water * x_replace - x_replace + L)
    A_1 = mu_oil / mu_water * A_2
    B_1 = p_center[time * 100 * t_step]
    B_2 = p_bound - A_2 * L

    x_1 = np.arange(0, x_replace, 0.001)
    x_2 = np.arange(x_replace, L, 0.001)

    for i in range(len(dist)):
        if dist[i] <= x_replace:
            p_oil = A_1 * dist[i] + B_1
            x_dist_pres[(set_coord[i], dist[i])].append(p_oil)
        else:
            p_water = A_2 * dist[i] + B_2
            x_dist_pres[(set_coord[i], dist[i])].append(p_water)

        x_dist_pres_exp[(set_coord[i], dist[i])].append(pt[i][time*100*t_step])

for key in x_dist_pres:
    plt.plot(x_dist_pres[key], 'r')
    plt.plot(x_dist_pres_exp[key])
    plt.title(key)
    plt.show()
