# рассчитывается давление в трещине по модели PKN с заданным давлением на скважине ( в центре трещины) и с миним. напряжением на кончике трещины.
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    l_fr0 = 0.2
    hx = 0.01
    hy = 0.01
    t_step = 1
    N_fr = int(l_fr0/hx)
    N_well = int(N_fr/2-1)+2
    nu = 0.2
    H = 0.07
    E = 3*10**9
    G = E/2/(1+nu)
    k = 4*(1-nu)*H/3.14/G
    mu = 0.1
    perm = 2*10**(-15)
    alpha = 1/24/mu/k/hx**2
    Pinj = 32*10**5
    Sh = 5*10**5
    Pres = 1*10**5
    w0 = k*(Pinj - Sh)
    coef = -perm/mu*(Pinj/2-Pres)/hy/2
    #q = np.ones((N_fr-1, 1))*coef
    q = np.zeros((N_fr-1, 1))

    T_exp = 100
    w = np.zeros((N_fr - 1, 1))

def Pressure_in_frac(N_fr, t_step, N_well, alpha, w0, q, w, k, Sh):

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

    for n in range(0, N_fr-1):
        if n+1 == N_well:
            A[n][n + 1] = 0
            B[n] = B[n] - alpha * (w[n] ** 3 + w[n+1]**3)*w0
        if n-1 == N_well:
            A[n][n - 1] = 0
            B[n] = B[n] - alpha * (w[n] ** 3 + w[n-1]**3)*w0

    A = np.delete(A, N_well, axis=0)
    A = np.delete(A, N_well, axis=1)
    B = np.delete(B, N_well)

    w_new = np.linalg.solve(A,B)
    w_new = np.insert(w_new, N_well, w0)


    w = w_new.reshape(N_fr-1, 1)
    P_new = w / k + Sh

    return P_new, w



if __name__ == '__main__':

    for t in range(T_exp):
        P_new, w_new = Pressure_in_frac(N_fr, t_step, N_well, alpha, w0, q, w, k, Sh)
        w = w_new

    print(P_new)
    fig = plt.figure()
    surf = plt.plot(P_new)
    plt.show()


