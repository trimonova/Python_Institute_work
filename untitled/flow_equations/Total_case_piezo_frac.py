import matplotlib.pyplot as plt
import numpy as np
from flow_in_frac import Pressure_in_frac
from matplotlib import cm
from piezo_for_wells_and_frac_new_BC import PorePressure_in_Time
from scipy import interpolate

from flow_equations.flow_in_frac_with_Q_in_well import Pressure_in_frac_with_Q_in_well

# входные данные для задачи течения в трещине
l_fr0 = 0.2     # начальные размер трещины
hx = 0.01    # шаг по х
hy = 0.01     # шаг по у
hz = 0.07
t_step = 0.01
N_fr_float = int(l_fr0/hx)
if N_fr_float%2 != 0:
    N_fr = N_fr_float+1
else:
    N_fr = N_fr_float
N_well = int(N_fr/2-1)
nu = 0.2
H = 0.07
E = 3*10**9
G = E/2/(1+nu)
k = 4*(1-nu)*H/3.14/G
mu = 0.001
perm = 2*10**(-15)
alpha_fr = 1/24/mu/k/hx**2
#Pinj = 13*10**5
Pinj = 0
Sh = 10*10**5
Pres = 1*10**5
w0 = k*(Pinj - Sh)

q = np.zeros((N_fr-1, 1))

Q = 0.0000001 # м3/с

w0_zasechki = 0.001 # раскрытие трещины у скважины
h_zasechki = 0.01
q[N_well] = Q/0.01/0.001 # нужно уточнять делители на эксперименте
print(q)

T_exp = 100

# входные данные для задачи пласта
alpha = 0.8*10**-12
beta = 0.17*10**-9

V = hx*hy*hz
coeff_1 = hx*hz/hy
coeff_2 = hy*hz/hx
Lx = 0.5
Ly = 0.5
N = int(Lx/hx) # количество ячеек вдоль оси х
M = int(Ly/hy)
if N%2 != 0:
    N = N+1
if M%2 != 0:
    M = M+1
print(N,M)

wells_with_Q = {}
wells_with_P = {(round(N/2+round(0.121/hx)),round(M/2+round(0.121/hy))): 15*10**5, (round(N/2-round(0.121/hx)),round(M/2-round(0.121/hy))): 1*10**5}
#frac_with_P = {(int(N/2), int(M/2)):25*10**5, (int(N/2)-1, int(M/2)):24*10**5, (int(N/2)+1, int(M/2)):24*10**5, (int(N/2)-2, int(M/2)):23*10**5, (int(N/2)+2, int(M/2)):23*10**5, (int(N/2)-3, int(M/2)):22*10**5, (int(N/2)+3, int(M/2)):22*10**5, (int(N/2)+4, int(M/2)):21*10**5, (int(N/2)-4, int(M/2)):21*10**5}
#wells_with_Q = {(int((Lx/2)/hx),int((Ly/2+0.121)/hy)): -0.00000003}
#wells_with_P = {(int((Lx/2+0.121)/hx),int((Ly/2+0.121)/hy)): 20*10**5, (int((Lx/2-0.121)/hx), int((Ly/2-0.121)/hy)): 1*10**5, (int((Lx/2-0.121)/hx), int((Ly/2)/hy)): 5*10**5, (int((Lx/2+0.121)/hx), int((Ly/2)/hy)): 5*10**5, (int((Lx/2)/hx), int((Ly/2-0.121)/hy)): 5*10**5, (int((Lx/2)/hx), int((Ly/2+0.121)/hy)): 5*10**5}
frac_with_P = {}
Pres = 1*10**5 # давление в пласте
w = np.zeros((N_fr - 1, 1))
Pres_distrib = np.ones((N, M)) * Pres

for t in range(T_exp):
    if Pinj != 0:
        P_in_frac, w_new =Pressure_in_frac(N_fr, t_step, N_well, alpha_fr, w0, q, w, k, Sh)
    else:
        P_in_frac, w_new = Pressure_in_frac_with_Q_in_well(N_fr, w, alpha_fr, t_step, q, k, Sh)

    for i in range(N_fr-1):
        frac_with_P[(int(N/2-N_fr/2+i+1), int(M/2))] = P_in_frac[i]

    P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, wells_with_P, frac_with_P, Pres, V, coeff_1, coeff_2, Pres_distrib)
    Pres_distrib = P_total
    w = w_new
    # считаем утечки
    Q_fr = []
    for i in range(N_fr-1):
        P_lost = P_total[int(N/2-N_fr/2+i+1)][int(M/2)+1]
        q_i = -perm/mu*(P_in_frac[i]-P_lost)/hy # возможно, нужно домножить на 2
        Q_fr.append(q_i)
    if Pinj == 0:
        Q_fr[N_well] = Q/0.01/0.001
        q = Q_fr
        q_sum = sum(q)
    else:
        Q_fr[N_well] = 0
        q = Q_fr
        q_sum = sum(q) - w0**2/12/mu*(P_total[N_well]-P_total[N_well-1])/hx
    print(q_sum)
print(w)


X = np.zeros((N,M))
Y = np.zeros((N, M))
for m in range(M):
    for n in range(N):
        X[n][m] = n*hx
        Y[n][m] = m*hy

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
P_list = [k for k in P_total.flat]


CP_list = zip(X_list, Y_list, P_list)

if __name__ == '__main__':
    print(min(P_list), max(P_list))

    xi = np.linspace(min(X_list),max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')

    levels = list(range(0,2700000,50000))
    fig = plt.figure()
    plt.subplot(211)
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)

    t = np.arange(0, 2 * np.pi, 0.01)
    r = 0.215
    plt.plot(r * np.sin(t) + Lx/2, r * np.cos(t) + Ly/2)

    #ax = fig.gca(projection='3d')

    #surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.subplot(212)
    plt.plot(P_in_frac)

    plt.show()

