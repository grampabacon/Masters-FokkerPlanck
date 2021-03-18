import numpy as np
import matplotlib.pyplot as plt
import math
import inputs.constants as const
import inputs.parameters as param

##############################################################
#                                                            #
# Functions used in the plotting file stationary_solution.py #
#                                                            #
##############################################################

# Range of mechanical energies for plotting
energy_min = 0
energy_max = 200
energy_divisions = 10000
mechanical_energies = np.linspace(energy_min, energy_max, energy_divisions + 1)
d_energy = (energy_max - energy_min) / energy_divisions

# Array of midpoint energies
midpoints = np.ndarray(len(mechanical_energies) - 1)
for j in range(len(mechanical_energies) - 1):
    midpoints[j] = (mechanical_energies[j] + mechanical_energies[j + 1]) / 2

# Theta for averaging over phase.
theta_divisions = 2000
thetas = np.linspace(0, np.pi, theta_divisions + 1)
d_theta = np.pi / theta_divisions

# Fermi distribution
def fermi_distribution(energy):
    return 1 / (1 + np.exp(energy / param.kT))


def gaussian(energy_change, epsilon_zero, width):
    return np.exp((energy_change - epsilon_zero) ** 2 / (2 * (width ** 2)))


def gaussian_norm():
    return 1


# Energy changes
def displacement(mechanical_energy, average=True):
    if average:
        return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energy * param.W_c) / param.oscillator_mass) * np.sin(thetas)
    else:
        return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energy * param.W_c) / param.oscillator_mass)


def energy_change_plus_left(mechanical_energy, average):
    return - param.W + param.W_L - (param.coupling_force * displacement(mechanical_energy, average) / param.W_c)


def energy_change_plus_right(mechanical_energy, average):
    return - param.W + param.W_R - (param.coupling_force * displacement(mechanical_energy, average) / param.W_c)


def energy_change_minus_left(mechanical_energy, average):
    return param.W - param.W_L + (param.coupling_force * displacement(mechanical_energy, average) / param.W_c)


def energy_change_minus_right(mechanical_energy, average):
    return param.W - param.W_R + (param.coupling_force * displacement(mechanical_energy, average) / param.W_c)


# Tunneling rates
def gamma_plus_right_modified(mechanical_energy, width, epsilon_zero, average=True):
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right(mechanical_energy, average)) * (1 - fermi_distribution(- energy_change_plus_right(mechanical_energy, average))) * (gaussian_norm() * gaussian(energy_change_minus_right(mechanical_energy, average), epsilon_zero, width))


def gamma_minus_right_modified(mechanical_energy, width, epsilon_zero, average=True):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(mechanical_energy, average)) * (fermi_distribution(energy_change_minus_right(mechanical_energy, average))) * (gaussian_norm() * gaussian(energy_change_minus_right(mechanical_energy, average), epsilon_zero, width))


def gamma_plus_left(mechanical_energy, average=True):
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left(mechanical_energy, average)) * (1 - fermi_distribution(- energy_change_plus_left(mechanical_energy, average)))


def gamma_minus_left(mechanical_energy, average=True):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(mechanical_energy, average)) * (fermi_distribution(energy_change_minus_left(mechanical_energy, average)))



# Derivatives of the tunneling rate with respect to W
def d_gamma_plus_right_modified(mechanical_energy, width, epsilon_zero, average=True):
    return ((2 * param.gamma_zero_right) / (width ** 2)) * np.exp((param.tunneling_exponent_right + (1 / param.kT)) * energy_change_minus_right(mechanical_energy, average)) * (gaussian(energy_change_minus_right(mechanical_energy, average), epsilon_zero, width) * gaussian_norm()) \
           * (- epsilon_zero + energy_change_minus_right(mechanical_energy, average) + (param.tunneling_exponent_right + (1 / param.kT)) * (width ** 2) + (- epsilon_zero + energy_change_minus_right(mechanical_energy, average) + param.tunneling_exponent_right * (width ** 2))) \
           * (fermi_distribution(energy_change_minus_right(mechanical_energy, average)) ** 2)


def d_gamma_minus_right_modified(mechanical_energy, width, epsilon_zero, average=True):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(mechanical_energy, average)) * (gaussian(energy_change_minus_right(mechanical_energy, average), epsilon_zero, width) * gaussian_norm()) * (- (1 / param.kT) * np.exp(energy_change_minus_right(mechanical_energy, average) / param.kT)
            + (1 + np.exp(energy_change_minus_right(mechanical_energy, average) / param.kT)) * (param.tunneling_exponent_right + ((- epsilon_zero + energy_change_minus_right(mechanical_energy, average)) / (width ** 2)))) \
            * (fermi_distribution(energy_change_minus_right(mechanical_energy, average)) ** 2)


def d_gamma_plus_left(mechanical_energy, average=True):
    return 2 * param.gamma_zero_left * np.exp((param.tunneling_exponent_left + (1 / param.kT)) * energy_change_minus_left(mechanical_energy)) * (param.tunneling_exponent_left + param.tunneling_exponent_left *
                                                                                                                                                 np.exp(energy_change_minus_left(mechanical_energy, average) / param.kT) +
                                                                                                                                                 (1 / param.kT)) * (fermi_distribution(energy_change_minus_left(mechanical_energy, average)) ** 2)

def d_gamma_minus_left(mechanical_energy, average=True):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(mechanical_energy)) * (param.tunneling_exponent_left + (param.tunneling_exponent_left - (1 / param.kT)) *
                                                                                                                          np.exp(energy_change_minus_left(mechanical_energy, average) / param.kT)) * (
                                                                                                                          fermi_distribution(energy_change_minus_left(mechanical_energy, average)) ** 2)


# Total on and off tunneling rates
def gamma_plus(mechanical_energy, width, epsilon_zero, average=True):
    return gamma_plus_left(mechanical_energy, average) + gamma_plus_right_modified(mechanical_energy, width, epsilon_zero, average)


def gamma_minus(mechanical_energy, width, epsilon_zero, average=True):
    return gamma_minus_left(mechanical_energy, average) + gamma_plus_right_modified(mechanical_energy, width, epsilon_zero, average)


# Total on and off tunneling rate derivatives
def gamma_plus_derivative_dw(mechanical_energy, width, epsilon_zero, average=True):
    return d_gamma_plus_left(mechanical_energy) + d_gamma_plus_right_modified(mechanical_energy, width, epsilon_zero, average)


def gamma_minus_derivative_dw(mechanical_energy, width, epsilon_zero, average=True):
    return d_gamma_minus_left(mechanical_energy) + d_gamma_minus_right_modified(mechanical_energy, width, epsilon_zero, average)


# Total tunneling rates
def gamma_total(mechanical_energy, width, epsilon_zero, average=True):
    return gamma_plus(mechanical_energy, width, epsilon_zero, average) + gamma_minus(mechanical_energy, width, epsilon_zero, average)


# Total tunneling rate derivative
def gamma_total_derivative_dw(mechanical_energy, width, epsilon_zero, average=True):
    return gamma_plus_derivative_dw(mechanical_energy, width, epsilon_zero, average) + gamma_minus_derivative_dw(mechanical_energy, width, epsilon_zero, average)


# Average island occupation
def average_occupation(mechanical_energy, width, epsilon_zero, average=True):
    return gamma_plus(mechanical_energy, width, epsilon_zero, average) / gamma_total(mechanical_energy, width, epsilon_zero, average)


# Occupation derivative
def average_occupation_derivative_dw(mechanical_energy, width, epsilon_zero, average=True):
    return (1 / gamma_total(mechanical_energy, width, epsilon_zero, average) ** 2) * ((gamma_minus(mechanical_energy, width, epsilon_zero, average) * gamma_plus_derivative_dw(mechanical_energy, width, epsilon_zero, average)) - (gamma_plus(mechanical_energy, width, epsilon_zero, average) * gamma_minus_derivative_dw(mechanical_energy, width, epsilon_zero, average)))


# # Kappa
# # Averaged over the oscillation phase
# def kappa1():
#     kappas = np.ndarray(len(mechanical_energies))
#     for i in range(len(mechanical_energies)):
#         dn = average_occupation_derivative_dw(mechanical_energies[i])
#         gamma_t = gamma_total(mechanical_energies[i])
#
#         _k = ((np.cos(thetas) ** 2) / np.pi) * (dn / gamma_t)
#         k_right = _k[1:]
#         k_left = _k[:-1]
#         k = (param.oscillator_frequency / param.quality_factor) + ((param.coupling_force * param.coupling_force / param.oscillator_mass)
#                                                                    * (d_theta / 2) * np.sum(k_left + k_right))
#
#         kappas[i] = k
#     return kappas


def kappa():
    kappas = np.ndarray(len(mechanical_energies))
    for i in range(len(mechanical_energies)):
        dn = average_occupation_derivative_dw(mechanical_energies[i])
        gamma_t = gamma_total(mechanical_energies[i])

        _k = ((np.cos(thetas) ** 2) / np.pi) * (dn / gamma_t)
        # k_right = _k[1:]
        # k_left = _k[:-1]
        # k = (param.oscillator_frequency / param.quality_factor) + ((param.coupling_force * param.coupling_force / param.oscillator_mass)
        #                                                            * (d_theta / 2) * np.sum(k_left + k_right))

        k = np.trapz(_k, dx=d_theta)
        k *= (param.coupling_force * param.coupling_force / param.oscillator_mass)
        k += param.oscillator_frequency / param.quality_factor

        kappas[i] = k
    return kappas


# Diffusion coefficient D(E)
# Averaged over the oscillation phase
def diffusion():
    diffusions = np.ndarray(len(mechanical_energies))
    for i in range(len(mechanical_energies)):
        n = average_occupation(mechanical_energies[i])
        gamma_t = gamma_total(mechanical_energies[i])

        _d = ((np.cos(thetas) ** 2) / np.pi) * (n * (1 - n)) / gamma_t
        # d_right = _d[1:]
        # d_left = _d[:-1]
        # d = (param.coupling_force ** 2 / param.oscillator_mass) * (d_theta / 2) * np.sum(d_left + d_right)
        d = np.trapz(_d, dx=d_theta)
        d *= (param.coupling_force ** 2 / param.oscillator_mass)

        diffusions[i] = d
    return diffusions
