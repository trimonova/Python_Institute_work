import numpy as np
import matplotlib.pyplot as plt

q = 0.000003
mu_oil = 0.11
mu_water = 0.003
perm = 2*10**(-15)
r_well = 0.0075
h = 0.07
por = 0.2
Q = q*mu_oil/perm/6.28/r_well/h
L = 0.215
Pb = 10**5
t = 10
a = q/6.28/r_well/h/por*t
deltaX_1 = a/2
deltaX_2 = (L-a)/2
P1 = Q*deltaX_1 + Pb + 2*Q*mu_water/mu_oil*deltaX_2
P2 = Pb + 2*Q*mu_water/mu_oil*deltaX_2

print(P1)
print(P2)
print(a)

plt.figure()
plt.scatter(0, P1)
plt.scatter(a, P2)
plt.scatter(L, Pb)

t = 20
a = q/6.28/r_well/h/por*t
deltaX_1 = a/2
deltaX_2 = (L-a)/2
P1 = Q*deltaX_1 + Pb + 2*Q*mu_water/mu_oil*deltaX_2
P2 = Pb + 2*Q*mu_water/mu_oil*deltaX_2

print(P1)
print(P2)
print(a)

plt.scatter(0, P1)
plt.scatter(a, P2)
plt.scatter(L, Pb)

plt.show()



