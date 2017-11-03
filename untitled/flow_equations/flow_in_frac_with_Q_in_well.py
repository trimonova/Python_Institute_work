# Рассчитывается давление в трещине при заданном расходе на скважину (в центре трещины), закачка задается в м/с, т.е возможно надо разделить закачку в м3/с на высоту
# трещины и на начальное раскрытие
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    l_fr0 = 0.2
    hx = 0.01
    hy = 0.01
    t_step = 0.01
    N_fr = int(l_fr0/hx)
    N_well = int(N_fr/2-1)
    nu = 0.2
    H = 0.07
    E = 3*10**9
    G = E/2/(1+nu)
    k = 4*(1-nu)*H/3.14/G
    mu = 0.001
    perm = 2*10**(-15)
    alpha = 1/24/mu/k/hx**2
    Pinj = 17*10**5
    Sh = 10*10**5
    Pres = 8*10**5
    w0 = k*(Pinj - Sh)
    coef = -perm/mu*(Pinj/2-Pres)/hy/1.5
    #q = np.ones((N_fr-1, 1))*coef
    q = np.zeros((N_fr-1, 1))

    w = np.zeros((N_fr-1, 1))
    q[N_well] = 0.001 # положительно значение - закачка
    print(q)
    T_exp = 500


def Pressure_in_frac_with_Q_in_well(N_fr, w, alpha, t_step, q, k, Sh):

    A = np.zeros((N_fr - 1, N_fr - 1))
    B = np.zeros((N_fr - 1, 1))
    for n in range(1, N_fr-2):

        A[n][n] = -alpha*(w[n+1]**3 + 2*w[n]**3 + w[n-1]**3) - 1/t_step
        A[n][n-1] = alpha*(w[n]**3 + w[n-1]**3)
        A[n][n+1] = alpha * (w[n] ** 3 + w[n+1]**3)

    A[0][0] = -alpha*(w[1]**3 + 2*w[0]**3) - 1/t_step
    A[0][1] = alpha*(w[0]**3 + w[1]**3)
    A[N_fr-2][N_fr-2] = -alpha*(w[N_fr-3]**3 + 2*w[N_fr-2]**3) - 1/t_step
    A[N_fr-2][N_fr-3] = alpha*(w[N_fr-2]**3 + w[N_fr-3]**3)

    for n in range(0, N_fr-1):
        B[n] = -1/t_step*w[n] - q[n]

    w_new = np.linalg.solve(A,B)

    w = w_new
    P_new = w / k + Sh

    return P_new, w_new

if __name__ == '__main__':
    for time in range(1,T_exp):
        P_new, w_new = Pressure_in_frac_with_Q_in_well(N_fr, w, alpha, t_step, q, k, Sh)
        w = w_new

    print(P_new)
    fig = plt.figure()
    surf = plt.plot(P_new)
    plt.show()


