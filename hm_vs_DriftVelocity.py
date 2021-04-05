import matplotlib.pyplot as plt
import numpy as np
from math import sqrt 

vm = range(0,0.1,10
close(mstring('all'))
clear(mstring('all'))
clc()
vm = mslice[0:0.1:10]
lm = 1
k = 1
w = 1
delta = 0.5
sigmam = delta / 2
lambdam = ((lm) ** 2) / sigmam


#part 1


l = 2. * (1 / sqrt(k)) *elmul* (k *elmul* vm - 1)
m = (-2) *elmul* (1 / sqrt(k)) *elmul* (k *elmul* vm + 1)
n = exp(8 * vm)
z1 = normcdf(l)
z2 = normcdf(m)

F1 = z1 + n *elmul* z2

k1 = k - 1
l1 = 2. * (1 / sqrt(k1)) *elmul* (k1 *elmul* vm - 1)
m1 = (-2) *elmul* (1 / sqrt(k1)) *elmul* (k1 *elmul* vm + 1)
n1 = exp(8 * vm)
z11 = normcdf(l1)
z21 = normcdf(m1)

F12 = z11 + n *elmul* z21

z3 = phi(1 / 0)
F2 = z3 + n *elmul* z3

F3 = F1 - F2
plot(vm, F3)
hold(mstring('on'))



for i in mslice[0:3]:
    k = k + 1
    l = 2. * (1 / sqrt(k)) *elmul* (k *elmul* vm - 1)
    m = (-2) *elmul* (1 / sqrt(k)) *elmul* (k *elmul* vm + 1)
    n = exp(8 * vm)
    z1 = normcdf(l)
    z2 = normcdf(m)

    F1 = z1 + n *elmul* z2

    k1 = k - 1
    l1 = 2. * (1 / sqrt(k1)) *elmul* (k1 *elmul* vm - 1)
    m1 = (-2) *elmul* (1 / sqrt(k1)) *elmul* (k1 *elmul* vm + 1)
    n1 = exp(8 * vm)
    z11 = normcdf(l1)
    z21 = normcdf(m1)

    F12 = z11 + n *elmul* z21

    z3 = phi(1 / 0)
    F2 = z3 + n *elmul* z3

    F31 = F1 - F12
    F3 = F1 - F2


    plot(vm, F31)
    hold(mstring('on'))

    legend(mstring('$h_m [n,1]$'), mstring('$h_m [n,2]$'), mstring('$h_m [n,3]$'), mstring('$h_m [n,4]$'), mstring('$h_m [n,5]$'), mstring('interpreter'), mstring('latex'))
    xlabel(mstring('$v_m[m/sec]$'), mstring('interpreter'), mstring('latex'))
    ylabel(mstring('$h_m [n,l]$'), mstring('interpreter'), mstring('latex'))

end



@mfunction("p")
def phi(z=None):
    p = 1 / 2. * (erfc(z / (sqrt(2))))
end