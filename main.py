import numpy as np
import matplotlib.pyplot as plt
import math
import inputs.constants as const
import inputs.parameters as param

# Theta for averaging over phase.
theta_divisions = 500
thetas = np.linspace(0, 2 * np.pi, theta_divisions)
d_theta = (2 * np.pi) / theta_divisions

# Range of mechanical energies for plotting
E_min = 0 * param.W_c
E_max = 1000 * param.W_c
energy_divisions = 10000
mechanical_energies = np.linspace(E_min, E_max, energy_divisions)
d_energy = (E_max - E_min) / energy_divisions

# Array of midpoint energies for averaging over theta
midpoints = np.ndarray(len(mechanical_energies) - 1)
for j in range(len(mechanical_energies) - 1):
    midpoints[j] = (mechanical_energies[j] + mechanical_energies[j + 1]) / 2


# Fermi distribution TODO: Check chemical potential
def fermi_distribution(energy):
    return 1 / (1 + np.exp(energy / (const.BOLTZMANN * param.temperature)))


def fermi_distribution_derivative_wrt_w_left(mechanical_energy):
    return (-np.exp(-energy_change_plus_left(mechanical_energy) / (const.BOLTZMANN * param.temperature)) / (const.BOLTZMANN * param.temperature)) * fermi_distribution(- energy_change_plus_left(mechanical_energy)) ** 2


def fermi_distribution_derivative_wrt_w_right(mechanical_energy):
    return (- (1 / (const.BOLTZMANN * param.temperature)) * np.exp(- energy_change_plus_right(mechanical_energy) / (const.BOLTZMANN * param.temperature))) * fermi_distribution(- energy_change_plus_right(mechanical_energy)) ** 2


# Energy changes
def displacement(mechanical_energy):
    return (1 / param.oscillator_frequency) * np.sqrt(np.sqrt((2 * mechanical_energy) / param.oscillator_mass)) * np.sin(thetas)


def energy_change_plus_left(mechanical_energy):
    return - param.W + param.W_L - param.coupling_force * displacement(mechanical_energy)


def energy_change_plus_right(mechanical_energy):
    return - param.W + param.W_R - param.coupling_force * displacement(mechanical_energy)


def energy_change_minus_left(mechanical_energy):
    return param.W - param.W_L + param.coupling_force * displacement(mechanical_energy)


def energy_change_minus_right(mechanical_energy):
    return param.W - param.W_R + param.coupling_force * displacement(mechanical_energy)


# Tunneling rates
def gamma_plus_left(mechanical_energy):
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left(mechanical_energy)) * (1 - fermi_distribution(- energy_change_plus_left(mechanical_energy)))


def gamma_plus_right(mechanical_energy):
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right(mechanical_energy)) * (1 - fermi_distribution(- energy_change_plus_right(mechanical_energy)))


def gamma_minus_left(mechanical_energy):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_plus_left(mechanical_energy)) * (fermi_distribution(energy_change_minus_left(mechanical_energy)))


def gamma_minus_right(mechanical_energy):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_plus_right(mechanical_energy)) * (fermi_distribution(energy_change_minus_right(mechanical_energy)))


# Tunneling rate derivatives
def gamma_plus_derivative_dw_left(mechanical_energy):
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left(mechanical_energy)) * (param.tunneling_exponent_left
                                                                                                                               - param.tunneling_exponent_left * fermi_distribution(- energy_change_plus_left(mechanical_energy))
                                                                                                                               - fermi_distribution_derivative_wrt_w_left(mechanical_energy))


def gamma_plus_derivative_dw_right(mechanical_energy):
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right(mechanical_energy)) * (param.tunneling_exponent_right
                                                                                                                                  - param.tunneling_exponent_right * fermi_distribution(- energy_change_plus_right(mechanical_energy))
                                                                                                                                  - fermi_distribution_derivative_wrt_w_right(mechanical_energy))


def gamma_minus_derivative_dw_left(mechanical_energy):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(mechanical_energy)) * (param.tunneling_exponent_left * fermi_distribution(energy_change_minus_left(mechanical_energy))
                                                                                                                          + fermi_distribution_derivative_wrt_w_left(mechanical_energy))


def gamma_minus_derivative_dw_right(mechanical_energy):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(mechanical_energy)) * (param.tunneling_exponent_right * fermi_distribution(energy_change_minus_right(mechanical_energy))
                                                                                                                             + fermi_distribution_derivative_wrt_w_right(mechanical_energy))


# Total on and off tunneling rates
def gamma_plus(mechanical_energy):
    return gamma_plus_left(mechanical_energy) + gamma_plus_right(mechanical_energy)


def gamma_minus(mechanical_energy):
    return gamma_minus_left(mechanical_energy) + gamma_minus_right(mechanical_energy)


# Total on and off tunneling rate derivatives
def gamma_plus_derivative_dw(mechanical_energy):
    return gamma_plus_derivative_dw_left(mechanical_energy) + gamma_plus_derivative_dw_right(mechanical_energy)


def gamma_minus_derivative_dw(mechanical_energy):
    return gamma_minus_derivative_dw_left(mechanical_energy) + gamma_minus_derivative_dw_right(mechanical_energy)


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
    return (1 / gamma_total(mechanical_energy) ** 2) * (gamma_total(mechanical_energy) * gamma_plus_derivative_dw(mechanical_energy) - gamma_plus(mechanical_energy) * gamma_total_derivative_dw(mechanical_energy))


# Kappa
# Averaged over the oscillation phase
def kappa():
    kappas = np.ndarray(len(midpoints))
    for i in range(len(midpoints)):
        dn = average_occupation_derivative_dw(midpoints[i])
        gamma_t = gamma_total(midpoints[i])

        _k = ((np.cos(thetas) ** 2) / np.pi) * (dn / gamma_t)
        k_right = _k[1:]
        k_left = _k[:-1]
        k = (param.oscillator_frequency / param.quality_factor) + ((param.coupling_force * param.coupling_force / param.oscillator_mass) * (d_theta / 2) * np.sum(k_left + k_right))

        kappas[i] = k
    return kappas


# Diffusion coefficient D(E)
# Averaged over the oscillation phase
def diffusion():
    diffusions = np.ndarray(len(midpoints))
    for i in range(len(midpoints)):
        n = average_occupation(midpoints[i])
        gamma_t = gamma_total(midpoints[i])

        _d = ((np.cos(thetas) ** 2) / np.pi) * (n * (1 - n)) / gamma_t
        d_right = _d[1:]
        d_left = _d[:-1]
        d = (param.coupling_force ** 2 / param.oscillator_mass) * midpoints[i] * (d_theta / 2) * np.sum(d_left + d_right)

        diffusions[i] = d
    return diffusions


# Calculating the integrals
def stationary_integrand():
    # return np.divide(np.multiply(-midpoints, kappa()), diffusion())
    return (-midpoints * kappa()) / diffusion()


def perform_stationary_integral(integrands):
    integrated = [0]
    for i in range(len(integrands) - 1):
        _integral = integrated[i]
        _integral += (integrands[i] + integrands[i + 1]) * (d_energy / 2)
        # math.floor(_integral)
        integrated.append(_integral)
    return integrated


############
# Plotting #
############
colors = ["r-", "b-", "g-", "m-", "c-", "y-", "k-", "r--", "b--", "g--", "m--"]


# Plot P(E) against W_c for  current parameters
def plot_p():
    plt.plot(midpoints / param.W_c, np.exp(perform_stationary_integral(stationary_integrand())), "k")
    plt.xlabel("Energy / W_c")
    plt.ylabel("P(E)")
    plt.show()


# Plots 4 regions on a single graph for different values of the bias voltage
def plots_regions_on_graph():
    plots = 4
    legend = [x for x in range(plots)]
    for i in range(plots):
        param.W = 5 * param.W_c
        param.bias_voltage = i * 2 * param.W_c / const.ELEMENTARY_CHARGE
        param.WL = - const.ELEMENTARY_CHARGE * param.bias_voltage / 2
        param.WR = const.ELEMENTARY_CHARGE * param.bias_voltage / 2

        legend[i] = "eV_b=" + str(param.bias_voltage * const.ELEMENTARY_CHARGE // param.W_c) + "W_c, W=" + str(param.W // param.W_c) + "W_c"
        plt.plot(midpoints / param.W_c, (np.exp(perform_stationary_integral(stationary_integrand()))), colors[i])
        print("Progress: " + str((i + 1) / plots * 100) + "%")
    plt.xlabel("Energy (W_c)")
    plt.ylabel("P(E)")
    plt.legend(legend)
    plt.show()


def plot_regions_seperate_graphs():
    plots = 4
    legend = [x for x in range(plots)]
    for i in range(plots):
        param.W = 0  # 5 * param.W_c
        param.bias_voltage = i * 2 * param.W_c / const.ELEMENTARY_CHARGE
        param.WL = - const.ELEMENTARY_CHARGE * param.bias_voltage / 2
        param.WR = const.ELEMENTARY_CHARGE * param.bias_voltage / 2

        plt.plot(midpoints / param.W_c, (np.exp(perform_stationary_integral(stationary_integrand()))), colors[i])
        plt.xlabel("Mechanical Energy / W_c")
        plt.ylabel("P(E)")
        plt.legend(["eV_b=" + str(param.bias_voltage * const.ELEMENTARY_CHARGE // param.W_c) + "W_c | W=" + str(param.W // param.W_c) + "W_c",])
        plt.show()
        print("Progress: " + str((i + 1) / plots * 100) + "%")


def plot_kappa():
    colors = ["r-", "b-", "g-", "m-", "c-", "y-", "k-", "r--", "b--", "g--", "m--"]
    for i in range(8):
        param.W = 0 * param.W_c
        param.bias_voltage = i * 2 * param.W_c / const.ELEMENTARY_CHARGE
        param.W_L = - const.ELEMENTARY_CHARGE * param.bias_voltage / 2
        param.W_R = const.ELEMENTARY_CHARGE * param.bias_voltage / 2
        plt.plot(midpoints / param.W_c, kappa(), colors[i])
    plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
    plt.xlabel("Energy [Wc]")
    plt.ylabel("Inverse damping time [Hz]")
    plt.legend(["eVb=2Wc, W=0", "eVb=4Wc, W=0", "eVb=6Wc, W=0", "eVb=8Wc, W=0"])
    plt.show()


# plot_regions_seperate_graphs()
plots_regions_on_graph()
# plot_p()
