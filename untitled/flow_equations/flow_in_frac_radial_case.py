# рассчитывается давление в трещине в радиальном случае в двух приближениях: большие утечки, маленькие утечки.
import numpy as np
if __name__ == '__main__':
    t = 2
    nu = 0.2
    mu = 0.001
    Q = 0.000001
    E = 3*10**9
    C = 0.0004
    S = 10*10**5
    G = E / 2 / (1 + nu)
    Rw = 0.0075 # для вертикальной трещины - длина интервала перфорации
    hx = 0.01
    print("hello")

def Pressure_in_frac_bigFluidLost(Q, t, C, nu, mu, E, G, S, Rw, hx):
    Ww = 1.92*((1-nu**2)**2*mu**2*Q**3/C/E**2)**(1/8)*t**(1/16)
    R = 1/3.14*Q**0.5/C**0.5*t**0.25
    f_rw = Rw / R
    Pw = S - 5/4/3.14*G*Ww/R*np.log10(f_rw)
    N_R = round(R/hx)
    P_in_frac = []
    for r in range(1,N_R+1):
        P_in_frac.append(Pw - 6*mu*Q/3.14/Ww**3*np.log10(r*hx/Rw))
    P_in_frac = np.hstack((list(reversed(P_in_frac)),Pw,P_in_frac))
    print(R,Pw)
    # print(P_in_frac)
    return P_in_frac, len(P_in_frac)

if __name__ == '__main__':
    Pressure_in_frac_bigFluidLost(Q, t, C, nu, mu, E, G, S, Rw, hx)

def Pressure_in_frac_NoFluidLost(Q, t, nu, mu, E, G, S, Rw, hx):
    Ww = 2.17*((1-nu**2)**2*mu**2*Q**3/E**2)**(1/9)*t**(1/9)
    R = 0.52*(E*Q**3/(1-nu**2)/mu)**(1/9)*t**(4/9)
    f_rw = Rw / R
    Pw = S - 5/4/3.14*G*Ww/R*np.log10(f_rw)
    N_R = round(R/hx)
    P_in_frac = []
    for r in range(1,N_R):
        P_in_frac.append(Pw - 6*mu*Q/3.14/Ww**3*np.log10(r*hx/Rw))
    P_in_frac = np.hstack((list(reversed(P_in_frac)),Pw,P_in_frac))
    print(R,Pw)
    print(P_in_frac)
    return P_in_frac, len(P_in_frac)

if __name__ == '__main__':
    P_in_frac, length = Pressure_in_frac_NoFluidLost(Q, t, nu, mu, E, G, S, Rw, hx)
    print(P_in_frac, length)

