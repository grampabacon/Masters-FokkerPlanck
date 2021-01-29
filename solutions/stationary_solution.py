import main
import numpy as np
import matplotlib.pyplot as plt
import math
import inputs.constants as const
import inputs.parameters as param

################################################
#                                              #
# See main.py for code used in these functions #
#                                              #
################################################


gate_capacitance = 0  # 10e-15
left_capacitance = 10e-15
right_capacitance = 10e-15


# Calculating the integrals
def stationary_integrand():
    # return (-main.mechanical_energies * main.kappa()) / main.diffusion()
    return - main.kappa() / main.diffusion()


def perform_stationary_integral(integrands):
    integrated = [0]
    for i in range(len(integrands) - 1):
        _i = integrated[i]
        _i += (integrands[i] + integrands[i + 1]) * (main.d_energy / 2)
        integrated.append(_i)
    return integrated


# def perform_stationary_integral1(integrands):
#     integrated = np.ndarray(len(integrands) - 1)
#     for i in range(len(integrands) - 1):
#         integrated[i] = (integrands[i] + integrands[i + 1]) * (main.d_energy / 2)
#     return integrated


############
# Plotting #
############
colors = ["r-", "b-", "g-", "m-", "c-", "y-", "k-", "r--", "b--", "g--", "m--"]


# Plot P(E) against W_c for  current parameters
# def plot_p():
#     plt.plot(main.mechanical_energies // param.W_c, np.exp(perform_stationary_integral(stationary_integrand())), "k")
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("P(E)")
#     plt.show()


# def plot_n():
#     # plt.plot(main.midpoints / param.W_c, (main.gamma_minus(main.midpoints) / main.gamma_total(main.midpoints)), "g")  # * np.exp(perform_stationary_integral(stationary_integrand()))
#     plt.plot(main.midpoints / param.W_c, (main.gamma_plus(main.midpoints) / main.gamma_total(main.midpoints)), "r")  # np.exp(perform_stationary_integral(stationary_integrand()))
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("n-bar")
#     plt.show()


# Plots 4 regions on a single graph for different values of the bias voltage
# def plots_regions_on_graph():
#     plots = 4
#     legend = [x for x in range(plots)]
#     for i in range(plots):
#         param.W = 5 * param.W_c
#         param.bias_voltage = i * 2 * param.W_c / const.ELEMENTARY_CHARGE
#         param.WL = - const.ELEMENTARY_CHARGE * param.bias_voltage / 2
#         param.WR = const.ELEMENTARY_CHARGE * param.bias_voltage / 2
#
#         legend[i] = "eV_b=" + str(param.bias_voltage * const.ELEMENTARY_CHARGE // param.W_c) + "W_c, W=" + str(param.W // param.W_c) + "W_c"
#         plt.plot(main.midpoints / param.W_c, (np.exp(perform_stationary_integral(stationary_integrand()))), colors[i])
#         print("Progress: " + str((i + 1) / plots * 100) + "%")
#     plt.xlabel("Energy (W_c)")
#     plt.ylabel("P(E)")
#     plt.legend(legend)
#     plt.show()


# def plot_regions_seperate_graphs():
#     plots = 1
#     for i in range(plots):
#         param.W = 1 * param.W_c  # 5 * param.W_c
#         param.bias_voltage = 1 * param.W_c / const.ELEMENTARY_CHARGE
#         param.WL = 0
#         param.WR = 2 * param.W_c
#
#         plt.plot(main.mechanical_energies / param.W_c, (np.exp(perform_stationary_integral(stationary_integrand()))), colors[i])
#         plt.xlabel("Mechanical Energy / W_c")
#         plt.ylabel("P(E)")
#         plt.legend(["eV_b=" + str(param.bias_voltage * const.ELEMENTARY_CHARGE // param.W_c) + "W_c | W=" + str(param.W // param.W_c) + "W_c", ])
#         plt.show()
#         print("Progress: " + str((i + 1) / plots * 100) + "%")


# def plot_multiple_kappa():
#     for i in range(4):
#         param.W = 0 * param.W_c
#         bias_voltage = i * 2 * param.W_c / const.ELEMENTARY_CHARGE
#         param.W_L = - const.ELEMENTARY_CHARGE * bias_voltage / 2
#         param.W_R = const.ELEMENTARY_CHARGE * bias_voltage / 2
#
#         plt.plot(main.midpoints / param.W_c, main.kappa(), colors[i])
#     plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("Kappa / Hz")
#     plt.legend(["eVb=2Wc, W=0", "eVb=4Wc, W=0", "eVb=6Wc, W=0", "eVb=8Wc, W=0"])
#     plt.show()


# def plot_diffusion():
#     param.W = 1 * param.W_c
#     param.W_L = 0  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
#     param.W_R = 2 * param.W_c
#
#     diffusion = main.diffusion()
#
#     plt.plot(main.mechanical_energies / param.W_c, diffusion, "k")
#     plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("Diffusion / J s^-1")
#     plt.title("W=" + str(param.W // param.W_c) + "W_c, Minimum diffusion=%.2f" % np.amin(diffusion))
#     plt.show()


# def plot_stationary_integrand():
#     param.W = 1 * param.W_c
#     param.W_L = 0  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
#     param.W_R = 2 * param.W_c
#
#     plt.plot(main.mechanical_energies / param.W_c, stationary_integrand(), "k")
#     plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("alpha")
#     plt.title("W=" + str(param.W // param.W_c) + "W_c")
#     plt.show()


# def plot_stationary_integral():
#     param.W = 1 * param.W_c
#     param.W_L = 0  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
#     param.W_R = 2 * param.W_c
#
#     plt.plot(main.midpoints / param.W_c, perform_stationary_integral(stationary_integrand()), "k")
#     plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("\gamma")
#     plt.title("W=" + str(param.W // param.W_c) + "W_c")
#     plt.show()


# def plot_dn_dw():
#     param.W = 1 * param.W_c
#     param.W_L = 0  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
#     param.W_R = 2 * param.W_c
#
#     dn_dw = np.ndarray(len(main.mechanical_energies))
#     for i in range(len(main.mechanical_energies)):
#         # dn_dw[i] = main.average_occupation_derivative_dw(main.mechanical_energies[i])
#         print(main.mechanical_energies[i])
#
#     plt.plot(main.mechanical_energies / param.W_c, dn_dw, "k")
#     plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("dn/dW")
#     plt.title("W=" + str(param.W // param.W_c) + "W_c, Minimum diffusion=%.2f" % np.amin(dn_dw))
#     plt.show()


# def plot_multiple_kappa():
#     legend = []
#     for i in range(2):
#         bias_voltage = i * 2 * param.W_c / const.ELEMENTARY_CHARGE
#         param.W_L = - const.ELEMENTARY_CHARGE * bias_voltage / 2
#         param.W_R = const.ELEMENTARY_CHARGE * bias_voltage / 2
#         for j in range(3):
#             param.W = j * param.W_c
#
#             legend.append("W=" + str(param.W // param.W_c) + ";eV_b=" + str(param.bias_voltage * const.ELEMENTARY_CHARGE // param.W_c))
#             plt.plot(main.mechanical_energies / param.W_c, main.kappa(), colors[i])
#     plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
#     plt.xlabel("Energy / W_c")
#     plt.ylabel("Kappa / Hz")
#     plt.legend(["eVb=2Wc, W=0", "eVb=4Wc, W=0", "eVb=6Wc, W=0", "eVb=8Wc, W=0"])
#     plt.show()


def normalise(probs):
    # s = np.linalg.norm(probs)
    s = np.sum(probs)
    print(s)
    normed = probs / s
    return normed


def plot_prob():
    param.W = 2
    bias = 9 / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
    param.W_R = w_right(bias)

    plt.plot(main.mechanical_energies, (np.exp(perform_stationary_integral(stationary_integrand()))), "k")
    plt.xlabel("Energy / W_c")
    plt.ylabel("P(E)")
    plt.title("$W$=" + str(round(param.W, 1)) + "$W_c$, $e V_b$=" + str(round(bias * const.ELEMENTARY_CHARGE, 1)) + "$W_c$, $W_L$=" + str(round(param.W_L, 1)) + "$W_c$, $W_R$=" + str(round(param.W_R, 1)) + "$W_c$", y=1.08)

    plt.show()


def plot_kappa():
    param.W = 2  # * param.W_c
    bias = 9 / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
    param.W_R = w_right(bias)

    kappas = main.kappa()

    plt.plot(main.mechanical_energies, kappas, "k")
    plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')
    plt.xlabel("Energy / W_c")
    plt.ylabel("Kappa / Hz")
    plt.title(
        "$W$=" + str(round(param.W, 1)) + "$W_c$, $e V_b$=" + str(round(bias * const.ELEMENTARY_CHARGE, 1)) + "$W_c$; $W_L$=" + str(round(param.W_L, 1)) + "$W_c$, $W_R$=" + str(round(param.W_R, 1)) + "$W_c$, Min($\kappa$)=%.2f" % np.amin(
            kappas), y=1.08)
    plt.show()


def find_roots():
    param.W = 2
    bias = 9 / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    kappas = main.kappa()
    arr = main.mechanical_energies * kappas

    signs = np.sign(arr)
    roots = []

    for i in range(len(signs)):
        if signs[i] == 0:
            roots.append((i * main.d_energy))
            continue
        if (i + 1) < len(signs):
            if signs[i + 1] != signs[i]:
                roots.append((i + (i + 1) / 2) * main.d_energy)
                continue
    return roots


def total_capacitance():
    return gate_capacitance + left_capacitance + right_capacitance


def w_left(bias):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / (2 * param.W_c)) - (right_capacitance * bias / const.ELEMENTARY_CHARGE))


def w_right(bias):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / (2 * param.W_c)) + (left_capacitance * bias / const.ELEMENTARY_CHARGE))


# plot_regions_seperate_graphs()
print(find_roots())
plot_prob()
