import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as colors
import inputs.constants as const
import inputs.parameters as param

mechanical_energy = 1000

theta_divisions = 300
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
    return (1 / total_capacitance()) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c)) - (right_capacitance * e_bias))


def w_right(e_bias):
    return (1 / total_capacitance()) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c)) + (left_capacitance * e_bias))


def e_bias_left(w):
    return (total_capacitance() / right_capacitance) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c * total_capacitance())) - w)


def e_bias_right(w):
    return (total_capacitance() / left_capacitance) * (w - (const.ELEMENTARY_CHARGE ** 2 / (2 * param.W_c * total_capacitance())))


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

    fig, ax = plt.subplots(num=f'E={mechanical_energy}Wc')

    # ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='RdGy', alpha=1)
    # im = ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='coolwarm', alpha=1)
    # ax.axis(aspect='image')

    divnorm = colors.TwoSlopeNorm(vcenter=0)

    im = ax.contourf(w_mesh, bias_mesh, kt, 20, cmap='plasma', norm=divnorm)  # Maybe plasma, viridis
    zero = ax.contour(w_mesh, bias_mesh, kt, 1, colors='black', levels=[0], linestyles='dashed')
    ax.clabel(zero, fmt='%1.1f')

    zero_contour = im.collections[0].get_paths()[0].vertices
    minimum = zero_contour[np.argmin(zero_contour, axis=0)[1]]

    label = "(" f'{float(f"{minimum[0]:.3g}"):g}' + ", " + f'{float(f"{minimum[1]:.3g}"):g}' + ")"
    ax.plot(minimum[0], minimum[1], "wp")
    # ax.annotate(label, (minimum[0], minimum[1]), textcoords="offset points", xytext=(1, -12), ha="center")

    handles = [mlines.Line2D([], [], color='w', marker='p', markersize=6, linewidth=0, markeredgecolor='black')]
    ax.legend(handles, [label], loc='lower right')

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(0, 10)

    ax.axhline(8, linestyle='solid', color='g', linewidth=0.3)
    # points = ax.plot([-2, -1.3, -0.5, 1], [4, 4, 4, 4], 'gs')
    point = ax.plot(-2.4, 8, 'gx')

    ax.tick_params(top=True, right=True)
    ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
    ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

    ax.set_title("Mechanical Energy: " + str(mechanical_energy) + "$W_c$")

    cb = fig.colorbar(im, orientation='vertical')
    cb.set_label("$\kappa$-tilde", size=14)
    # cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    ax.set_aspect(0.75)

    filename = "KappaDensity " + f"E={mechanical_energy}Wc. Min+Mark" + ".png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    fig.show()


def plot_kappa_tilde_symmetric():
    w_list2 = np.linspace(-5, 5, 200)
    e_bias_list2 = np.linspace(-10, 10, 200)
    w_mesh2, bias_mesh2 = np.meshgrid(w_list2, e_bias_list2)

    kt2 = kappa_tilde_mesh(w_mesh2, bias_mesh2)

    fig, ax = plt.subplots()

    # ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='RdGy', alpha=1)
    # im = ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='coolwarm', alpha=1)
    # ax.axis(aspect='image')

    divnorm = colors.TwoSlopeNorm(vcenter=0)

    im = ax.contourf(w_mesh2, bias_mesh2, kt2, 20, cmap='plasma', norm=divnorm) # Maybe plasma, viridis
    ax.contour(w_mesh2, bias_mesh2, kt2, 1, colors='black', levels=[0], linestyles='dashed')

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(-10, 10)

    ax.tick_params(top=True, right=True)
    ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
    ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

    ax.set_title("Mechanical Energy: " + str(mechanical_energy) + "$W_c$")

    cb = fig.colorbar(im, orientation='vertical')
    cb.set_label("$\kappa$-tilde", size=14)
    # cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    ax.set_aspect(0.4)

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


def plot_kappa_tildes():
    global mechanical_energy

    x, y = -0.5, 2.2
    for i in range(4):
        mechanical_energy = (i + 0) * 20

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

        point = ax.plot(x, y, 'mx')
        handles = [mlines.Line2D([], [], color='m', marker='x', markersize=6, linewidth=0)]
        ax.legend(handles, ["(" + str(x) + ", " + str(y) + ")"], loc='lower right')

        ax.tick_params(top=True, right=True)
        ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
        ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

        ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
        ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

        ax.set_title("Mechanical Energy: " + str(mechanical_energy) + "$W_c$")

        cb = fig.colorbar(im, orientation='vertical')
        cb.set_label("$\kappa$-tilde", size=14)
        # cb.ax.set_title("$\partial_W n$", size=14)
        cb.ax.tick_params(labelsize='large')

        ax.set_aspect(0.75)
        fig.show()


def plot_kappa_tilde_region_iii():
    global mechanical_energy

    mechanical_energy = 0
    kt0 = kappa_tilde_mesh(w_mesh, bias_mesh)

    _kt = np.sign(kt0)

    sign_changes = np.zeros((len(w_list), len(e_bias_list)))
    for i in range(250):
        mechanical_energy = (i + 1) * 20

        kt = np.sign(kappa_tilde_mesh(w_mesh, bias_mesh))
        for j in range(len(w_list)):
            for k in range(len(e_bias_list)):
                if _kt[j][k] != kt[j][k]:
                    sign_changes[j][k] += 1
        _kt = kt
    for j in range(len(w_list)):
        for k in range(len(e_bias_list)):
            if sign_changes[j][k] == 0:
                continue
            if sign_changes[j][k] == 1:
                sign_changes[j][k] = 0

    blue_colors = [(0, 0, 1, c) for c in np.linspace(0, 1, 100)]
    cmapblue = colors.LinearSegmentedColormap.from_list('mycmap', blue_colors, N=(np.amax(sign_changes) + 1))

    fig, ax = plt.subplots(num='Region iii')

    divnorm = colors.TwoSlopeNorm(vcenter=0)

    im = ax.contourf(w_mesh, bias_mesh, kt0, 20, cmap='plasma', norm=divnorm)  # Maybe plasma, viridis
    ax.contour(w_mesh, bias_mesh, kt0, 1, colors='black', levels=[0], linestyles='dashed')

    regions = ax.contourf(w_mesh, bias_mesh, sign_changes, 20, cmap=cmapblue)

    ax.axhline(3, linestyle='solid', color='g', linewidth=0.3)
    point = ax.plot([-2, -1.3, -0.5, 1], [3, 3, 3, 3], 'gs')

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(0, 10)

    # handles = [mlines.Line2D([], [], color='m', marker='x', markersize=6, linewidth=0)]
    # ax.legend(handles, ["(" + str(x) + ", " + str(y) + ")"], loc='lower right')

    ax.tick_params(top=True, right=True)
    ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
    ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

    # ax.set_title("Mechanical Energy: " + str(mechanical_energy) + "$W_c$")

    cb_regions = fig.colorbar(regions, ticks=np.linspace(0, 5, 6))
    cb_regions.set_label("Roots")
    cb_regions.ax.tick_params(labelsize='large')

    cb = fig.colorbar(im, orientation='vertical')
    cb.set_label("$\kappa$-tilde", size=14)
    # cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    ax.set_aspect(0.75)
    fig.show()


plot_kappa_tilde()