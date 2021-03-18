import solutions.modified_main as mmain
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import math
import inputs.constants as const
import inputs.parameters as param
mpl.use('Qt5Agg')

################################################
#                                              #
# See main.py for code used in these functions #
#                                              #
################################################


gate_capacitance = 0  # 10e-15
left_capacitance = 10e-15
right_capacitance = 10e-15


# Calculating the integrals
def stationary_integrand(kappas):
    # return (-main.mechanical_energies * main.kappa()) / main.diffusion()
    return - kappas / mmain.diffusion()


def perform_stationary_integral(integrands):
    integrated = [0]
    for i in range(len(integrands) - 1):
        _i = integrated[i]
        _i += (integrands[i] + integrands[i + 1]) * (mmain.d_energy / 2)
        integrated.append(_i)
    return integrated


############
# Plotting #
############
colors = ["r-", "b-", "g-", "m-", "c-", "y-", "k-", "r--", "b--", "g--", "m--"]


def normalise(probs):
    s = np.linalg.norm(probs)
    # s = np.sum(probs)
    # s = np.trapz(probs, main.mechanical_energies)
    print("Normalisation constant: " + str(s))
    normed = probs / s
    return normed


def plot_prob(bias, kappas):
    # param.W = 2
    # bias = 9 / const.ELEMENTARY_CHARGE
    # param.W_L = w_left(bias)  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
    # param.W_R = w_right(bias)

    plot1 = plt.figure("Probability")
    plt.plot(mmain.mechanical_energies, normalise(np.exp(perform_stationary_integral(stationary_integrand(kappas)))), "k")

    plt.xlabel("Energy / W_c")
    plt.ylabel("P(E)")
    plt.title("$W$=" + str(round(param.W, 2)) + "$W_c$, $e V_b$=" + str(round(bias * const.ELEMENTARY_CHARGE, 2)) + "$W_c$, $W_L$=" + str(round(param.W_L, 2)) + "$W_c$, $W_R$=" + str(round(param.W_R, 2)) + "$W_c$", y=1.03)

    roots = find_roots(kappas)

    handles = [mpl_patches.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="white", linewidth=0, alpha=0)]
    labels = ["Extrema at: [" + ", ".join(map(str, roots)) + "]$W_c$"]
    plt.legend(handles, labels, loc='best', fontsize='large', fancybox=True, framealpha=0.7, handlelength=0, handletextpad=0)

    filename = "Probability W=" + str(round(param.W, 2)) + ", eVb=" + str(round(bias * const.ELEMENTARY_CHARGE, 2)) + ", WL=" + str(round(param.W_L, 2)) + ", WR=" + str(round(param.W_R, 2)) + ".png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    plt.show()


def plot_kappa(bias, kappas):
    # param.W = 2  # * param.W_c
    # bias = 9 / const.ELEMENTARY_CHARGE
    # param.W_L = w_left(bias)  # - const.ELEMENTARY_CHARGE * bias_voltage / 2
    # param.W_R = w_right(bias)

    plot2 = plt.figure("Kappa")
    plt.plot(mmain.mechanical_energies, kappas, "k")
    plt.axhline(y=0, color='k', linewidth=0.5, label='_nolegend_')

    plt.xlabel("Energy / W_c")
    plt.ylabel("Kappa / Hz")
    plt.title(
        "$W$=" + str(round(param.W, 2)) + "$W_c$, $e V_b$=" + str(round(bias * const.ELEMENTARY_CHARGE, 2)) + "$W_c$; $W_L$=" + str(round(param.W_L, 2)) + "$W_c$, $W_R$=" + str(round(param.W_R, 2)), y=1.03)

    roots = find_roots(kappas)
    # for i in range(len(roots)):
    #     plt.axvline(roots[i], color='r')

    handles = [mpl_patches.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="white", linewidth=0, alpha=0)]
    labels = []
    if len(roots) == 0:
        labels.append("No roots.")
    else:
        labels.append("Roots: [" + ", ".join(map(str, roots)) + "]$W_c$")
    plt.legend(handles, labels, loc='best', fontsize='large', fancybox=True, framealpha=0.7, handlelength=0, handletextpad=0)

    filename = "Kappa W=" + str(round(param.W, 2)) + ", eVb=" + str(round(bias * const.ELEMENTARY_CHARGE, 2)) + ", WL=" + str(round(param.W_L, 1)) + ", WR=" + str(round(param.W_R, 2)) + ".png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    plt.show()


def plot_rates_plus_minus(w, e_bias):
    param.W = w
    bias = e_bias / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    width = 100
    epsilon_zero = 1

    fig, ax = plt.subplots(num='Plus Minus Rates')

    gamma_plus = mmain.gamma_plus(mmain.mechanical_energies, width, epsilon_zero, False)
    gamma_minus = mmain.gamma_minus(mmain.mechanical_energies, width, epsilon_zero, False)

    ax.plot(mmain.mechanical_energies, gamma_plus, 'g', label='$\Gamma^+$', linestyle='--')
    ax.plot(mmain.mechanical_energies, gamma_minus, 'r', label='$\Gamma^-$', linestyle='-.')

    ax.xaxis.set_label_text("E (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    ax.tick_params(axis='x', direction='in', labelsize=12)
    ax.tick_params(axis='y', direction='in', labelsize=12)

    ax.legend()

    title = "Total Tunneling Rates onto and off the Island"
    ax.set_title(title)

    fig.show()


def plot_rates_plus_left_right_subplots(w, e_bias):
    param.W = w
    bias = e_bias / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    width = 1
    epsilon_zero = 1

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), num='Plus Rates Modified Subplots')

    gamma_left = mmain.gamma_plus_left(mmain.mechanical_energies, False)
    gamma_right = mmain.gamma_plus_right_modified(mmain.mechanical_energies, width, epsilon_zero, False)

    axes[0].plot(mmain.mechanical_energies, gamma_left, 'g', label='$\Gamma^+_L$', linestyle='--')
    axes[1].plot(mmain.mechanical_energies, gamma_right, 'r', label='$\Gamma^+_R$', linestyle='-.')

    axes[0].xaxis.set_label_text("E (units of $W_c$)", fontsize=14)
    axes[1].xaxis.set_label_text("E (units of $W_c$)", fontsize=14)

    axes[0].yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    axes[0].tick_params(axis='x', direction='in', labelsize=12)
    axes[0].tick_params(axis='y', direction='in', labelsize=12)

    axes[1].tick_params(axis='x', direction='in', labelsize=12)
    axes[1].tick_params(axis='y', direction='in', labelsize=12)

    axes[0].legend()
    axes[1].legend()

    params = "W=" + str(w) + "$W_c$, $e V_b$=" + str(e_bias) + "$W_c$, $\epsilon_0$=" + str(epsilon_zero) + "$W_c$, $\sigma$=" + str(width) + "$W_c$"
    title = "Modified Tunneling Rates onto the Island\n" + params
    fig.suptitle(title)

    fig.tight_layout()


def plot_rates_plus_left_right_single(w, e_bias):
    param.W = w
    bias = e_bias / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    width = 1
    epsilon_zero = 1

    fig, ax = plt.subplots(num='Plus Rates Modified Single')

    gamma_left = mmain.gamma_plus_left(mmain.mechanical_energies, False)
    gamma_right = mmain.gamma_plus_right_modified(mmain.mechanical_energies, width, epsilon_zero, False)

    ax.plot(mmain.mechanical_energies, gamma_left, 'g', label='$\Gamma^+_L$', linestyle='--')
    ax.plot(mmain.mechanical_energies, gamma_right, 'r', label='$\Gamma^+_R$', linestyle='-.')

    ax.xaxis.set_label_text("E (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    ax.legend()

    params = "W=" + str(w) + "$W_c$, $e V_b$=" + str(e_bias) + "$W_c$, $\epsilon_0$=" + str(epsilon_zero) + "$W_c$, $\sigma$=" + str(width) + "$W_c$"
    title = "Modified Tunneling Rates off of the Island\n" + params
    ax.set_title(title)

    fig.tight_layout()
    
    
def plot_rates_minus_left_right_subplots(w, e_bias):
    param.W = w
    bias = e_bias / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    width = 1
    epsilon_zero = 1

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), num='minus Rates Modified Subplots')

    gamma_left = mmain.gamma_minus_left(mmain.mechanical_energies, False)
    gamma_right = mmain.gamma_minus_right_modified(mmain.mechanical_energies, width, epsilon_zero, False)

    axes[0].plot(mmain.mechanical_energies, gamma_left, 'g', label='$\Gamma^+_L$', linestyle='--')
    axes[1].plot(mmain.mechanical_energies, gamma_right, 'r', label='$\Gamma^+_R$', linestyle='-.')

    axes[0].xaxis.set_label_text("E (units of $W_c$)", fontsize=14)
    axes[1].xaxis.set_label_text("E (units of $W_c$)", fontsize=14)

    axes[0].yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    axes[0].tick_params(axis='x', direction='in', labelsize=12)
    axes[0].tick_params(axis='y', direction='in', labelsize=12)

    axes[1].tick_params(axis='x', direction='in', labelsize=12)
    axes[1].tick_params(axis='y', direction='in', labelsize=12)

    axes[0].legend()
    axes[1].legend()

    params = "W=" + str(w) + "$W_c$, $e V_b$=" + str(e_bias) + "$W_c$, $\epsilon_0$=" + str(epsilon_zero) + "$W_c$, $\sigma$=" + str(width) + "$W_c$"
    title = "Modified Tunneling Rates onto the Island\n" + params
    fig.suptitle(title)

    fig.tight_layout()


def plot_rates_minus_left_right_single(w, e_bias):
    param.W = w
    bias = e_bias / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    width = 1
    epsilon_zero = 1

    fig, ax = plt.subplots(num='Minus Rates Modified Single')

    gamma_left = mmain.gamma_minus_left(mmain.mechanical_energies, False)
    gamma_right = mmain.gamma_minus_right_modified(mmain.mechanical_energies, width, epsilon_zero, False)

    ax.plot(mmain.mechanical_energies, gamma_left, 'g', label='$\Gamma^-_L$', linestyle='--')
    ax.plot(mmain.mechanical_energies, gamma_right, 'r', label='$\Gamma^-_R$', linestyle='-.')

    ax.xaxis.set_label_text("E (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    ax.legend()

    params = "W=" + str(w) + "$W_c$, $e V_b$=" + str(e_bias) + "$W_c$, $\epsilon_0$=" + str(epsilon_zero) + "$W_c$, $\sigma$=" + str(width) + "$W_c$"
    title = "Modified Tunneling Rates off of the Island\n" + params
    ax.set_title(title)

    fig.show()


def find_roots(kappas):
    # param.W = 2
    # bias = 9 / const.ELEMENTARY_CHARGE
    # param.W_L = w_left(bias)
    # param.W_R = w_right(bias)

    signs = np.sign(kappas)
    roots = []

    for i in range(len(signs)):
        if signs[i] == 0:
            roots.append((i * mmain.d_energy))
            continue
        if (i + 1) < len(signs):
            if signs[i + 1] != signs[i]:
                roots.append(i * mmain.d_energy)
                continue
    if not roots:
        roots.append(0)
    return np.around(roots, 1)


def total_capacitance():
    return gate_capacitance + left_capacitance + right_capacitance


def w_left(bias):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / (2 * param.W_c)) - (right_capacitance * bias / const.ELEMENTARY_CHARGE))


def w_right(bias):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / (2 * param.W_c)) + (left_capacitance * bias / const.ELEMENTARY_CHARGE))


def plot_kappa_single(w, e_bias):
    param.W = w
    bias = e_bias / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    kappas = mmain.kappa()

    plot_kappa(bias, kappas)


# Plot the main graphs and data related to P(E). This plot will print the roots where the maxima of the probability function are located and the kappa(E) of the output P(E).
def plot_graphs(w, e_bias):
    param.W = w
    bias = e_bias / const.ELEMENTARY_CHARGE
    param.W_L = w_left(bias)
    param.W_R = w_right(bias)

    kappas = mmain.kappa()

    plot_kappa(bias, kappas)
    plot_prob(bias, kappas)


plot_rates_minus_left_right_subplots(0, -5)
plot_rates_minus_left_right_single(0, -5)
# plot_rates_plus_left_right_subplots(0, -5)
# plot_rates_plus_left_right_single(0, -5)
