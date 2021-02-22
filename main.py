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
energy_max = 50000
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


# def fermi_distribution_derivative_wrt_w_left(mechanical_energy):
#     return (- np.exp(- energy_change_plus_left(mechanical_energy) / param.kT) / param.kT) * (fermi_distribution(- energy_change_plus_left(mechanical_energy)) ** 2)
#
#
# def fermi_distribution_derivative_wrt_w_right(mechanical_energy):
#     return (- np.exp(- energy_change_plus_right(mechanical_energy) / param.kT) / param.kT) * (fermi_distribution(- energy_change_plus_right(mechanical_energy)) ** 2)


# Energy changes
def displacement(mechanical_energy):
    return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energy * param.W_c) / param.oscillator_mass) * np.sin(thetas)


def energy_change_plus_left(mechanical_energy):
    return - param.W + param.W_L - (param.coupling_force * displacement(mechanical_energy) / param.W_c)


def energy_change_plus_right(mechanical_energy):
    return - param.W + param.W_R - (param.coupling_force * displacement(mechanical_energy) / param.W_c)


def energy_change_minus_left(mechanical_energy):
    return param.W - param.W_L + (param.coupling_force * displacement(mechanical_energy) / param.W_c)


def energy_change_minus_right(mechanical_energy):
    return param.W - param.W_R + (param.coupling_force * displacement(mechanical_energy) / param.W_c)


# Tunneling rates
def gamma_plus_left(mechanical_energy):
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left(mechanical_energy)) * (1 - fermi_distribution(- energy_change_plus_left(mechanical_energy)))


def gamma_plus_right(mechanical_energy):
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right(mechanical_energy)) * (1 - fermi_distribution(- energy_change_plus_right(mechanical_energy)))


def gamma_minus_left(mechanical_energy):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(mechanical_energy)) * (fermi_distribution(energy_change_minus_left(mechanical_energy)))


def gamma_minus_right(mechanical_energy):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(mechanical_energy)) * (fermi_distribution(energy_change_minus_right(mechanical_energy)))


# # Tunneling rate derivatives
# def gamma_plus_derivative_dw_left(mechanical_energy):
#     return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left(mechanical_energy)) * (param.tunneling_exponent_left
#                                                                                                                                - param.tunneling_exponent_left * fermi_distribution(- energy_change_plus_left(mechanical_energy))
#                                                                                                                                - fermi_distribution_derivative_wrt_w_left(mechanical_energy))
#
#
# def gamma_plus_derivative_dw_right(mechanical_energy):
#     return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right(mechanical_energy)) * (param.tunneling_exponent_right
#                                                                                                                                   - param.tunneling_exponent_right * fermi_distribution(- energy_change_plus_right(mechanical_energy))
#                                                                                                                                   - fermi_distribution_derivative_wrt_w_right(mechanical_energy))
#
#
# def gamma_minus_derivative_dw_left(mechanical_energy):
#     return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(mechanical_energy)) * (param.tunneling_exponent_left * fermi_distribution(energy_change_minus_left(mechanical_energy))
#                                                                                                                           + fermi_distribution_derivative_wrt_w_left(mechanical_energy))
#
#
# def gamma_minus_derivative_dw_right(mechanical_energy):
#     return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(mechanical_energy)) * (param.tunneling_exponent_right * fermi_distribution(energy_change_minus_right(mechanical_energy))
#                                                                                                                              + fermi_distribution_derivative_wrt_w_right(mechanical_energy))


# Derivatives of the tunneling rate with respect to W
def d_gamma_plus_left(mechanical_energy):
    return 2 * param.gamma_zero_left * np.exp((param.tunneling_exponent_left + (1 / param.kT)) * energy_change_minus_left(mechanical_energy)) * (param.tunneling_exponent_left + param.tunneling_exponent_left *
                                                                                                                                                 np.exp(energy_change_minus_left(mechanical_energy) / param.kT) +
                                                                                                                                                 (1 / param.kT)) * (fermi_distribution(energy_change_minus_left(mechanical_energy)) ** 2)


def d_gamma_plus_right(mechanical_energy):
    return 2 * param.gamma_zero_right * np.exp((param.tunneling_exponent_right + (1 / param.kT)) * energy_change_minus_right(mechanical_energy)) * (param.tunneling_exponent_right + param.tunneling_exponent_right *
                                                                                                                                                    np.exp(energy_change_minus_right(mechanical_energy) / param.kT) +
                                                                                                                                                    (1 / param.kT)) * (fermi_distribution(energy_change_minus_right(mechanical_energy)) ** 2)


def d_gamma_minus_left(mechanical_energy):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(mechanical_energy)) * (param.tunneling_exponent_left + (param.tunneling_exponent_left - (1 / param.kT)) *
                                                                                                                          np.exp(energy_change_minus_left(mechanical_energy) / param.kT)) * (
                                                                                                                          fermi_distribution(energy_change_minus_left(mechanical_energy)) ** 2)


def d_gamma_minus_right(mechanical_energy):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(mechanical_energy)) * (param.tunneling_exponent_right + (param.tunneling_exponent_right - (1 / param.kT)) *
                                                                                                                             np.exp(energy_change_minus_right(mechanical_energy) / param.kT)) * (
                                                                                                                             fermi_distribution(energy_change_minus_right(mechanical_energy)) ** 2)


# Total on and off tunneling rates
def gamma_plus(mechanical_energy):
    return gamma_plus_left(mechanical_energy) + gamma_plus_right(mechanical_energy)


def gamma_minus(mechanical_energy):
    return gamma_minus_left(mechanical_energy) + gamma_minus_right(mechanical_energy)


# Total on and off tunneling rate derivatives
def gamma_plus_derivative_dw(mechanical_energy):
    return d_gamma_plus_left(mechanical_energy) + d_gamma_plus_right(mechanical_energy)


def gamma_minus_derivative_dw(mechanical_energy):
    return d_gamma_minus_left(mechanical_energy) + d_gamma_minus_right(mechanical_energy)


# Total tunneling rates
def gamma_total(mechanical_energy):
    return gamma_plus(mechanical_energy) + gamma_minus(mechanical_energy)


# Total tunneling rate derivative
def gamma_total_derivative_dw(mechanical_energy):
    return gamma_plus_derivative_dw(mechanical_energy) + gamma_minus_derivative_dw(mechanical_energy)


# Average island occupation
def average_occupation(mechanical_energy):
    return gamma_plus(mechanical_energy) / gamma_total(mechanical_energy)


# Occupation derivative
def average_occupation_derivative_dw(mechanical_energy):
    return (1 / gamma_total(mechanical_energy) ** 2) * ((gamma_minus(mechanical_energy) * gamma_plus_derivative_dw(mechanical_energy)) - (gamma_plus(mechanical_energy) * gamma_minus_derivative_dw(mechanical_energy)))


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
