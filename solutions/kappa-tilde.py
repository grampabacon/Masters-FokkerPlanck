import numpy as np
import matplotlib.pyplot as plt
import inputs.constants as const
import inputs.parameters as param

mechanical_energy = 0

theta_divisions = 200
theta_list = np.linspace(0, np.pi, theta_divisions)
d_theta = np.pi / theta_divisions

w_list = np.linspace(-5, 5, 200)
e_bias_list = np.linspace(0, 10, 200)
w_mesh, bias_mesh = np.meshgrid(w_list, e_bias_list)

gate_capacitance = 0  # 10e-15
left_capacitance = 10e-15
right_capacitance = 10e-15


# Fermi-Dirac Distribution with \mu = 0
def fermi_distribution(energy_change):
    return 1 / (1 + np.exp(energy_change / param.kT))


# Calculates oscillator transverse displacement using the mechanical energy of the oscillator
def displacement():
    return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energy * param.W_c) / param.oscillator_mass) * np.sin(theta_list)


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
    return (1 / total_capacitance()) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c)) - (right_capacitance * e_bias)) - (param.coupling_force * displacement() / param.W_c)


def w_right(e_bias):
    return (1 / total_capacitance()) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c)) + (left_capacitance * e_bias)) - (param.coupling_force * displacement() / param.W_c)


def e_bias_left(w):
    return (total_capacitance() / right_capacitance) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c * total_capacitance())) - w - (param.coupling_force * displacement() / param.W_c))


def e_bias_right(w):
    return (total_capacitance() / left_capacitance) * (w + (param.coupling_force * displacement() / param.W_c) - (const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c * total_capacitance())))


def kappa_tilde_mesh(w, e_bias):
    kappa_tilde_m = np.ndarray((len(w_list), len(e_bias_list)))

    for i in range(len(w_list)):
        for j in range(len(e_bias_list)):
            dn = d_n(w[i][j], e_bias[i][j])
            gamma_t = gamma_total(w[i][j], e_bias[i][j])

            _term = dn / gamma_t
            _k = ((np.cos(theta_list) ** 2) / np.pi) * (dn / gamma_t)

            k = np.trapz(_k, dx=d_theta)
            k *= (param.coupling_force * param.coupling_force / param.oscillator_mass)

            kappa_tilde_m[i][j] = k
    return kappa_tilde_m


def kappa_tilde_bias_slice(w, e_bias):
    kappa_tilde_b = np.ndarray(len(w_list))

    for i in range(len(w_list)):
        dn = d_n(w[i], e_bias)
        gamma_t = gamma_total(w[i], e_bias)

        _term = dn / gamma_t
        _k = ((np.cos(theta_list) ** 2) / np.pi) * (dn / gamma_t)

        k = np.trapz(_k, dx=d_theta)
        k *= (param.coupling_force * param.coupling_force / param.oscillator_mass)

        kappa_tilde_b[i] = k
    return kappa_tilde_b


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
    cb.set_label("$\partial_W n$", size=14)
    # cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    ax.set_aspect(0.75)
    fig.show()


def plot_kappa_tilde():
    kt = kappa_tilde_mesh(w_mesh, bias_mesh)

    fig, ax = plt.subplots()

    # ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='RdGy', alpha=1)
    # im = ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='coolwarm', alpha=1)
    # ax.axis(aspect='image')

    im = ax.contourf(w_mesh, bias_mesh, kt, 20, cmap='coolwarm')
    ax.contour(w_mesh, bias_mesh, kt, 1, colors='black', levels=[0], linestyles='dashed')

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(0, 10)

    ax.tick_params(top=True, right=True)
    ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
    ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

    cb = fig.colorbar(im, orientation='vertical')
    cb.set_label("$\kappa$-tilde", size=14)
    # cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    ax.set_aspect(0.75)
    fig.show()


def plot_kappa_tilde_dn_dw_bias_slice():
    kt = kappa_tilde_bias_slice(w_list, 6)
    #kt = kt / np.linalg.norm(kt)
    dn = d_n(w_list, 6)
    #dn = dn / np.linalg.norm(dn)

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('$W$ (units of $W_c$)', fontsize=14)
    ax1.set_ylabel('$\kappa$-tilde ($s^{-1}$)', color=color, fontsize=14)
    ax1.plot(w_list, kt, color=color)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('$\partial_W n$ ($W_c^{-1}$)', color=color, fontsize=14)
    ax2.plot(w_list, dn, color=color)
    ax2.tick_params(axis='y', labelcolor=color, labelsize=12)

    ax1.set_title("Comparing $\partial_W n$ and $\kappa$-tilde for $e V_b = 6 W_c$", fontsize=14)
    plt.tight_layout()
    fig.show()


plot_kappa_tilde_dn_dw_bias_slice()
