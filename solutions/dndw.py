import numpy as np
import matplotlib.pyplot as plt
import inputs.constants as const
import inputs.parameters as param

mechanical_energy = 1000

w_list = np.linspace(-5, 5, 2000)
e_bias_list = np.linspace(0, 10, 2000)
w_mesh, bias_mesh = np.meshgrid(w_list, e_bias_list)

gate_capacitance = 0  # 10e-15
left_capacitance = 10e-15
right_capacitance = 10e-15


# Fermi-Dirac Distribution with \mu = 0
def fermi_distribution(energy_change):
    return 1 / (1 + np.exp(energy_change / param.kT))


# Calculates oscillator transverse displacement using the mechanical energy of the oscillator
def displacement():
    return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energy * param.W_c) / param.oscillator_mass)


# Changes in energy available for tunneling after an electron tunnels on or off the island from the left or the right
def energy_change_plus_left(w, e_bias):
    return - w + w_left(e_bias) - (param.coupling_force * displacement() / param.W_c)


def energy_change_plus_right(w, e_bias):
    return - w + w_right(e_bias) - (param.coupling_force * displacement() / param.W_c)


def energy_change_minus_left(w, e_bias):
    return w - w_left(e_bias) + (param.coupling_force * displacement() / param.W_c)


def energy_change_minus_right(w, e_bias):
    return w - w_right(e_bias) + (param.coupling_force * displacement() / param.W_c)


# Tunneling rates
def gamma_plus_left(w, e_bias):
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left(w, e_bias)) * (1 - fermi_distribution(- energy_change_plus_left(w, e_bias)))


def gamma_plus_right(w, e_bias):
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right(w, e_bias)) * (1 - fermi_distribution(- energy_change_plus_right(w, e_bias)))


def gamma_minus_left(w, e_bias):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(w, e_bias)) * (fermi_distribution(energy_change_minus_left(w, e_bias)))


def gamma_minus_right(w, e_bias):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(w, e_bias)) * (fermi_distribution(energy_change_minus_right(w, e_bias)))


# Derivatives of the tunneling rate with respect to W
def d_gamma_plus_left(w, e_bias):
    return 2 * param.gamma_zero_left * np.exp((param.tunneling_exponent_left + (1 / param.kT)) * energy_change_minus_left(w, e_bias)) * (param.tunneling_exponent_left +
                                                                                                                                         param.tunneling_exponent_left * np.exp(energy_change_minus_left(w, e_bias) / param.kT) +
                                                                                                                                         (1 / param.kT)) * (fermi_distribution(energy_change_minus_left(w, e_bias)) ** 2)


def d_gamma_plus_right(w, e_bias):
    return 2 * param.gamma_zero_right * np.exp((param.tunneling_exponent_right + (1 / param.kT)) * energy_change_minus_right(w, e_bias)) * (param.tunneling_exponent_right +
                                                                                                                                            param.tunneling_exponent_right * np.exp(energy_change_minus_right(w, e_bias) / param.kT) +
                                                                                                                                            (1 / param.kT)) * (fermi_distribution(energy_change_minus_right(w, e_bias)) ** 2)


def d_gamma_minus_left(w, e_bias):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(w, e_bias)) * (param.tunneling_exponent_left +
                                                                                                                  (param.tunneling_exponent_left - (1 / param.kT)) *
                                                                                                                  np.exp(energy_change_minus_left(w, e_bias) / param.kT)) * (fermi_distribution(energy_change_minus_left(w, e_bias)) ** 2)


def d_gamma_minus_right(w, e_bias):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(w, e_bias)) * (param.tunneling_exponent_right +
                                                                                                                     (param.tunneling_exponent_right - (1 / param.kT)) *
                                                                                                                     np.exp(energy_change_minus_right(w, e_bias) / param.kT)) * (fermi_distribution(energy_change_minus_right(w, e_bias)) ** 2)


def gamma_plus(w, e_bias):
    return gamma_plus_left(w, e_bias) + gamma_plus_right(w, e_bias)


def gamma_minus(w, e_bias):
    return gamma_minus_left(w, e_bias) + gamma_minus_right(w, e_bias)


def d_gamma_plus(w, e_bias):
    return d_gamma_plus_left(w, e_bias) + d_gamma_plus_right(w, e_bias)


def d_gamma_minus(w, e_bias):
    return d_gamma_minus_left(w, e_bias) + d_gamma_minus_right(w, e_bias)


def gamma_total(w, e_bias):
    return gamma_plus(w, e_bias) + gamma_minus(w, e_bias)


def d_n(w, e_bias):
    return (1 / (gamma_total(w, e_bias) ** 2)) * (gamma_minus(w, e_bias) * d_gamma_plus(w, e_bias) - gamma_plus(w, e_bias) * d_gamma_minus(w, e_bias))


def total_capacitance():
    return gate_capacitance + left_capacitance + right_capacitance


def w_left(e_bias):
    return (1 / total_capacitance()) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c)) - (right_capacitance * e_bias))


def w_right(e_bias):
    return (1 / total_capacitance()) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c)) + (left_capacitance * e_bias))


def e_bias_left(w):
    return (total_capacitance() / right_capacitance) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c * total_capacitance())) - w)


def e_bias_right(w):
    return (total_capacitance() / left_capacitance) * (w - (const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c * total_capacitance())))


# def kappa_tilde(w, e_bias):
#     kappa_tilde_mesh = np.ndarray((len(w_list), len(e_bias_list)))
#
#     dn = d_n(w, e_bias, True)
#     print(dn)
#     gamma_t = gamma_total(w, e_bias, True)
#
#     for i in range(len(w_list)):
#         for j in range(len(e_bias_list)):
#             _term = dn[i][j] / gamma_t[i][j]
#             _k = ((np.cos(theta_mesh) ** 2) / np.pi) * (dn / gamma_t)
#
#             k = np.trapz(_k, dx=d_theta)
#             k *= (param.coupling_force * param.coupling_force / param.oscillator_mass)
#
#             kappa_tilde_mesh[i][j] = k
#     return kappa_tilde_mesh


def plot_dn_dw():
    dn = d_n(w_mesh, bias_mesh)

    fig, ax = plt.subplots()

    # ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='RdGy', alpha=1)
    # im = ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='coolwarm', alpha=1)
    # ax.axis(aspect='image')

    im = ax.contourf(w_mesh, bias_mesh, dn, 20, cmap='coolwarm')
    ax.contour(w_mesh, bias_mesh, dn, 1, colors='black', levels=[0], linestyles='dashed')

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(0, 10)

    ax.tick_params(top=True, right=True)
    ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
    ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

    cb = fig.colorbar(im, orientation='vertical')
    # cb.set_label("$\partial_W n$", size=14)
    cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    ax.set_aspect(0.75)
    fig.show()


plot_dn_dw()
