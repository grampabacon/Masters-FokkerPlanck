import numpy as np
import matplotlib.pyplot as plt
import inputs.parameters as param
import inputs.constants as const

gate_capacitance = 10e-15
left_capacitance = 10e-15
right_capacitance = 10e-15


# Range of mechanical energies for plotting
energy_min = 0 * param.W_c
energy_max = 1000 * param.W_c
energy_divisions = 10000
mechanical_energies = np.linspace(energy_min, energy_max, energy_divisions + 1)
d_energy = (energy_max - energy_min) / energy_divisions

# Array of midpoint energies
midpoints = np.ndarray(len(mechanical_energies) - 1)
for j in range(len(mechanical_energies) - 1):
    midpoints[j] = (mechanical_energies[j] + mechanical_energies[j + 1]) / 2


# Fermi distribution TODO: Check chemical potential
def fermi_distribution(energy):
    return 1 / (1 + np.exp(energy / (const.BOLTZMANN * param.temperature)))


def fermi_distribution_derivative_wrt_w_left():
    return (- np.exp(- energy_change_plus_left() / (const.BOLTZMANN * param.temperature)) / (const.BOLTZMANN * param.temperature)) * fermi_distribution(- energy_change_plus_left()) ** 2


def fermi_distribution_derivative_wrt_w_right():
    return (- np.exp(- energy_change_plus_right() / (const.BOLTZMANN * param.temperature)) / (const.BOLTZMANN * param.temperature)) * fermi_distribution(- energy_change_plus_right()) ** 2


# Energy changes
def displacement():
    return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energies) / param.oscillator_mass)


def energy_change_plus_left():
    return - param.W + param.W_L - param.coupling_force * displacement()


def energy_change_plus_right():
    return - param.W + param.W_R - param.coupling_force * displacement()


def energy_change_minus_left():
    return param.W - param.W_L + param.coupling_force * displacement()


def energy_change_minus_right():
    return param.W - param.W_R + param.coupling_force * displacement()


# Tunneling rates
def gamma_plus_left():
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left()) * (1 - fermi_distribution(- energy_change_plus_left()))


def gamma_plus_right():
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right()) * (1 - fermi_distribution(- energy_change_plus_right()))


def gamma_minus_left():
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_plus_left()) * (fermi_distribution(energy_change_minus_left()))


def gamma_minus_right():
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_plus_right()) * (fermi_distribution(energy_change_minus_right()))


# Tunneling rate derivatives
def gamma_plus_derivative_dw_left():
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left()) * (param.tunneling_exponent_left
                                                                                                              - param.tunneling_exponent_left * fermi_distribution(- energy_change_plus_left())
                                                                                                              - fermi_distribution_derivative_wrt_w_left())


def gamma_plus_derivative_dw_right():
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right()) * (param.tunneling_exponent_right
                                                                                                                 - param.tunneling_exponent_right * fermi_distribution(- energy_change_plus_right())
                                                                                                                 - fermi_distribution_derivative_wrt_w_right())


def gamma_minus_derivative_dw_left():
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left()) * (param.tunneling_exponent_left * fermi_distribution(energy_change_minus_left())
                                                                                                         + fermi_distribution_derivative_wrt_w_left())


def gamma_minus_derivative_dw_right():
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right()) * (param.tunneling_exponent_right * fermi_distribution(energy_change_minus_right())
                                                                                                            + fermi_distribution_derivative_wrt_w_right())


# Total on and off tunneling rates
def gamma_plus():
    return gamma_plus_left() + gamma_plus_right()


def gamma_minus():
    return gamma_minus_left() + gamma_minus_right()


# Total on and off tunneling rate derivatives
def gamma_plus_derivative_dw():
    return gamma_plus_derivative_dw_left() + gamma_plus_derivative_dw_right()


def gamma_minus_derivative_dw():
    return gamma_minus_derivative_dw_left() + gamma_minus_derivative_dw_right()


# Total tunneling rates
def gamma_total():
    return gamma_plus() + gamma_minus()


# Total tunneling rate derivative
def gamma_total_derivative_dw():
    return gamma_plus_derivative_dw() + gamma_minus_derivative_dw()


# Average island occupation
def average_occupation():
    return gamma_plus() / gamma_total()


# Occupation derivative
def average_occupation_derivative_dw():
    return (1 / (gamma_total() ** 2)) * ((gamma_minus() * gamma_plus_derivative_dw()) - (gamma_plus() * gamma_minus_derivative_dw()))


def total_capacitance():
    return gate_capacitance + left_capacitance + right_capacitance


def w_left(bias):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / 2) - (right_capacitance * bias / const.ELEMENTARY_CHARGE))


def w_right(bias):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / 2) + (left_capacitance * bias / const.ELEMENTARY_CHARGE))


def check_n_dn_dw():
    param.W = param.W_c

    bias = 0 * param.W_c / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    plt.plot(mechanical_energies / param.W_c, average_occupation_derivative_dw())

    n = average_occupation()
    print(n)
    plotted = n * (1 - n) / (const.BOLTZMANN * param.temperature)

    plt.plot(mechanical_energies / param.W_c, plotted * param.W_c, "g-")
    plt.legend(["dn/dW", "$n(1-n)/k_BT$"])
    plt.title("$V_b=0$; $W=$" + str(round(param.W / param.W_c, 1)) + "$W_c$")
    plt.xlabel("Energy / $W_c$")
    plt.ylabel("$n$")
    plt.show()


def plot_n():
    param.W = 0 * param.W_c

    bias = 15 * param.W_c / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    plt.plot(mechanical_energies / param.W_c, average_occupation())
    plt.title("$eV_b=$" + str(round(bias * const.ELEMENTARY_CHARGE / param.W_c, 1)) + "$W_c$; $W=$" + str(round(param.W / param.W_c, 1)) + "$W_c$; $W_L$=" + str(round(param.W_L / param.W_c, 1)) + "$W_c$, $W_R$=" + str(round(param.W_R / param.W_c, 1)) + "$W_c$", y=1.08)
    plt.xlabel("Energy / Wc")
    plt.ylabel("Average Island Occupation")
    plt.show()


plot_n()
