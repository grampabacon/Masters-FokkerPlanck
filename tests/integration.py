import numpy as np
from scipy import integrate

x_min = 0
x_max = 90
x_divisions = 1000000
x_list = np.linspace(x_min, x_max, x_divisions + 1)
d_x = (x_max - x_min) / x_divisions


def f(x):
    return x ** 2


def perform_integral(integrands):
    integrated = [0]
    for i in range(len(integrands) - 1):
        _i = integrated[i]
        _i += (integrands[i] + integrands[i + 1]) * (d_x / 2)
        integrated.append(_i)
    return integrated


def p(x):
    integrated = np.ndarray(len(x) - 1)
    for i in range(len(x) - 1):
        integrated[i] = (x[i] + x[i + 1]) * (d_x / 2)
    return integrated


def perform_cumtrap(integrands):
    return integrate.cumtrapz(x_list, integrands, dx=d_x)


def perform_trap(integrands):
    return np.trapz(integrands, dx=d_x)


#print(np.sum(perform_cumtrap((x_list))))
print(np.sum(p(f(x_list))))
# print(perform_trap(f(x_list)))
