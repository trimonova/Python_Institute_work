import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from piezo_for_wells_and_frac_new_BC import PorePressure_in_Time
from scipy import interpolate

from flow_equations.flow_in_frac_radial_case import Pressure_in_frac_bigFluidLost

# входные данные для задачи течения в трещине
l_fr0 = 0.2     # начальные размер трещины
hx = 0.01    # шаг по х
hy = 0.01     # шаг по у
hz = 0.07
t_step = 0.1
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

Sh = 10*10**5
Pres = 1*10**5

Q = 0.000001 # м3/с
C = 0.0001
T_exp = 300
Rw = 0.0075 # для вертикальной трещины - длина интервала перфорации
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

Pres_distrib = np.ones((N, M)) * Pres

for t in range(1,T_exp):

    P_in_frac, N_fr =Pressure_in_frac_bigFluidLost(Q, t*t_step, C, nu, mu, E, G, Sh, Rw, hx)

    N_well = int((N_fr-1)/2)
    print(P_in_frac, N_fr, N_well)

    for i in range(N_fr):
        frac_with_P[(int(N/2-N_fr/2+i+1), int(M/2))] = P_in_frac[i]

    P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, wells_with_P, frac_with_P, Pres, V, coeff_1, coeff_2, Pres_distrib)
    Pres_distrib = P_total

    # считаем утечки
    Q_fr = []
    for i in range(N_fr):
        P_lost = P_total[int(N / 2 - N_fr / 2 + i + 1)][int(M / 2) + 1]
        q_i = -perm / mu * (P_in_frac[i] - P_lost) / hy  # возможно, нужно домножить на 2
        Q_fr.append(q_i)

    #print(Q_fr)
    q = Q_fr
    try:
        Q_fr[N_well] = Q/0.001/0.01
    except IndexError:
        Q_fr.append(Q)

    q_sum = sum(q)

    print(q_sum)



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

