# Здесь рассматривается прямолинейное поршневое вытеснение маслом воды (Извеков + Басниев Кочина, стр. 204).
#За Нестационарность вычислений отвечает меняющиеся граничное условие на скважине. Но здесь нестационарность рассматривается только косвенно.
# Давление берется постоянным в интересующий нас момент времени. Можем брать несколько моментов времени, в них будет отличаться давления на забое, можно смотреть, как будут отличаться графики.
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
Texp = 2000 #
t_step = 1

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
for time in range(1,Texp+1):

    delta_p = p_bound - p_center[time*100*t_step]
    D = ((L*mu_water - (mu_oil - mu_water)*x_replace)/(mu_oil - mu_water))**2 + 4*(L*mu_water*x_replace - perm/fi*delta_p*t_step)/(mu_oil - mu_water)
    a_1 = (-(L*mu_water - (mu_oil - mu_water)*x_replace)/(mu_oil - mu_water) + D**0.5)/2
    a_2 = (-(L*mu_water - (mu_oil - mu_water)*x_replace)/(mu_oil - mu_water) - D**0.5)/2
    print(a_2, a_1)
    x_replace = max(a_2, a_1)
    print(p_center[time*100*t_step])

print(x_replace)

print(set_with_P)

A_2 = delta_p/(mu_oil/mu_water*x_replace - x_replace + L)
A_1 = mu_oil/mu_water * A_2
B_1 = p_center[Texp*100*t_step]
B_2 = p_bound - A_2*L

x_1 = np.arange(0, x_replace, 0.0001)
x_2 = np.arange(x_replace, L, 0.0001)

p_oil = A_1*x_1 + B_1
p_water = A_2*x_2 + B_2

fig = plt.figure()
plt.plot(x_1, p_oil)
plt.plot(x_2, p_water)

for i in range(np.shape(y)[1]):
    plt.scatter(dist[i], set_with_P[(x[0][i], y[0][i])][Texp*100*t_step])
    print(dist[i], set_with_P[(x[0][i], y[0][i])][Texp*100*t_step])
    i += 1
plt.show()