import main
import numpy as np
import matplotlib.pyplot as plt
import math
import inputs.constants as const
import inputs.parameters as param

# Theta for averaging over phase.
theta_divisions = 1000
thetas = np.linspace(0, np.pi, theta_divisions + 1)
d_theta = (np.pi) / theta_divisions


def fermi_distribution(energy):
    return 1 / (1 + np.exp(energy / (const.BOLTZMANN * 8)))


def fermi_distribution_derivatives(energy):
    # return
    derivatives = np.ndarray(len(energy))
    for i in range(len(energy)):
        print((- np.exp(energy[i] / (const.BOLTZMANN * 8)) / (const.BOLTZMANN * 8)) * fermi_distribution(energy[i]) ** 2)
        derivatives[i] = (- np.exp(energy[i] / (const.BOLTZMANN * 8)) / (const.BOLTZMANN * 8)) * fermi_distribution(energy[i]) ** 2
    return derivatives


def fermi_distribution_derivative(energy):
    return (- np.exp(energy / (const.BOLTZMANN * 8)) / (const.BOLTZMANN * 8)) * fermi_distribution(energy) ** 2


def plot_fermi():
    energies = np.linspace(-12500 * param.W_c, 12500 * param.W_c)
    plt.plot(energies, fermi_distribution(energies))
    plt.show()


def plot_fermi_derivative():
    energies = np.linspace(1, 125000 * param.W_c, 10000)
    print(fermi_distribution_derivatives(energies))
    plt.plot(energies, fermi_distribution_derivatives(energies))
    plt.show()


print(fermi_distribution_derivative(10000000 * 5E-3 * 1.6E-19))
