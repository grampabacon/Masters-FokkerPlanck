import main
import numpy as np
import matplotlib.pyplot as plt
import math
import inputs.constants as const
import inputs.parameters as param


# Calculating the integrals
def stationary_integrand():
    # return np.divide(np.multiply(-midpoints, kappa()), diffusion())
    return (-main.midpoints * main.kappa()) / main.diffusion()


def perform_stationary_integral(integrands):
    integrated = [0]
    for i in range(len(integrands) - 1):
        _integral = integrated[i]
        _integral += (integrands[i] + integrands[i + 1]) * (main.d_energy / 2)
        # math.floor(_integral)
        integrated.append(_integral)
    return integrated


############
# Plotting #
############
colors = ["r-", "b-", "g-", "m-", "c-", "y-", "k-", "r--", "b--", "g--", "m--"]


# Plot P(E) against W_c for  current parameters
def plot_p():
    plt.plot(main.midpoints / param.W_c, np.exp(perform_stationary_integral(stationary_integrand())), "k")
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
        plt.plot(main.midpoints / param.W_c, (np.exp(perform_stationary_integral(stationary_integrand()))), colors[i])
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

        plt.plot(main.midpoints / param.W_c, (np.exp(perform_stationary_integral(stationary_integrand()))), colors[i])
        plt.xlabel("Mechanical Energy / W_c")
        plt.ylabel("P(E)")
        plt.legend(["eV_b=" + str(param.bias_voltage * const.ELEMENTARY_CHARGE // param.W_c) + "W_c | W=" + str(param.W // param.W_c) + "W_c", ])
        plt.show()
        print("Progress: " + str((i + 1) / plots * 100) + "%")


def plot_kappa():
    for i in range(3):
        param.W = 0 * param.W_c
        param.bias_voltage = i * param.W_c / const.ELEMENTARY_CHARGE
        param.W_L = - const.ELEMENTARY_CHARGE * param.bias_voltage / 2
        param.W_R = const.ELEMENTARY_CHARGE * param.bias_voltage / 2
        plt.plot(main.midpoints / param.W_c, main.kappa(), colors[i])
    plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
    plt.xlabel("Energy / W_c")
    plt.ylabel("Kappa / Hz")
    plt.legend(["eVb=2Wc, W=0", "eVb=4Wc, W=0", "eVb=6Wc, W=0", "eVb=8Wc, W=0"])
    plt.show()


# plot_regions_seperate_graphs()
# plots_regions_on_graph()
plot_p()
