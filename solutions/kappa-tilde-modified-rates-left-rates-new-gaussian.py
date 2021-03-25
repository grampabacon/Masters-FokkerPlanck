import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as colors
import inputs.constants as const
import inputs.parameters as param
mpl.rcParams['figure.dpi'] = 150


mechanical_energy = 0

theta_divisions = 400
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


def gaussian(energy_change, epsilon_zero, width):
    return np.exp(- (energy_change - epsilon_zero) ** 2 / (2 * (width ** 2)))


def norm_constant():
    return 1


# def gamma_plus_left_e0(epsilon_zero):
#     return 1  # 2 * param.gamma_zero_left * np.exp(param.tunneling_exponent_left * epsilon_zero) * (1 - fermi_distribution(epsilon_zero))
#
#
# def gamma_minus_left_e0(epsilon_zero):
#     return 1  # param.gamma_zero_left * np.exp(param.tunneling_exponent_left * epsilon_zero) * (fermi_distribution(epsilon_zero))


# Calculates oscillator transverse displacement using the mechanical energy of the oscillator
def displacement(average):
    if average:
        return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energy * param.W_c) / param.oscillator_mass) * np.sin(theta_list)
    else:
        return (1 / param.oscillator_frequency) * np.sqrt((2 * mechanical_energy * param.W_c) / param.oscillator_mass)


# Changes in energy available for tunneling after an electron tunnels on or off the island from the left or the right
def energy_change_plus_left(w, e_bias, average):
    return - w + w_left(e_bias) - (param.coupling_force * displacement(average) / param.W_c)


def energy_change_plus_right(w, e_bias, average):
    return - w + w_right(e_bias) - (param.coupling_force * displacement(average) / param.W_c)


def energy_change_minus_left(w, e_bias, average):
    return w - w_left(e_bias) + (param.coupling_force * displacement(average) / param.W_c)


def energy_change_minus_right(w, e_bias, average):
    return w - w_right(e_bias) + (param.coupling_force * displacement(average) / param.W_c)


# Tunneling rates
def gamma_plus_left_modified(w, e_bias, epsilon_zero, width, average=True):
    return 2 * param.gamma_zero_left * np.exp(- param.tunneling_exponent_left * energy_change_plus_left(w, e_bias, average)) * (1 - fermi_distribution(- energy_change_plus_left(w, e_bias, average))) * (norm_constant() * gaussian(energy_change_minus_left(w, e_bias, average), epsilon_zero, width))


def gamma_minus_left_modified(w, e_bias, epsilon_zero, width, average=True):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(w, e_bias, average)) * (fermi_distribution(energy_change_minus_left(w, e_bias, average))) * (norm_constant() * gaussian(energy_change_minus_left(w, e_bias, average), epsilon_zero, width))


def gamma_plus_right(w, e_bias, average=True):
    return 2 * param.gamma_zero_right * np.exp(- param.tunneling_exponent_right * energy_change_plus_right(w, e_bias, average)) * (1 - fermi_distribution(- energy_change_plus_right(w, e_bias, average)))


def gamma_minus_right(w, e_bias, average=True):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(w, e_bias, average)) * (fermi_distribution(energy_change_minus_right(w, e_bias, average)))


# Derivatives of the tunneling rate with respect to W
def d_gamma_plus_left_modified(w, e_bias, epsilon_zero, width, average=True):
    return ((2 * param.gamma_zero_left) / (width ** 2)) * np.exp((param.tunneling_exponent_left + (1 / param.kT)) * energy_change_minus_left(w, e_bias, average)) * (norm_constant() * gaussian(energy_change_minus_left(w, e_bias, average), epsilon_zero, width)) \
           * ((1 / param.kT) * (width ** 2) + (1 + np.exp(energy_change_minus_left(w, e_bias, average))) * (epsilon_zero - energy_change_minus_left(w, e_bias, average) + param.tunneling_exponent_left * (width ** 2))) \
           * (fermi_distribution(energy_change_minus_left(w, e_bias, average)) ** 2)


def d_gamma_minus_left_modified(w, e_bias, epsilon_zero, width, average=True):
    return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(w, e_bias, average)) * (norm_constant() * gaussian(energy_change_minus_left(w, e_bias, average), epsilon_zero, width)) * (- (1 / param.kT) * np.exp(energy_change_minus_left(w, e_bias, average) / param.kT)
            + (1 + np.exp(energy_change_minus_left(w, e_bias, average) / param.kT)) * (param.tunneling_exponent_left - ((- epsilon_zero + energy_change_minus_left(w, e_bias, average)) / (width ** 2)))) \
            * (fermi_distribution(energy_change_minus_left(w, e_bias, average)) ** 2)


def d_gamma_plus_right(w, e_bias, average=True):
    return 2 * param.gamma_zero_right * np.exp((param.tunneling_exponent_right + (1 / param.kT)) * energy_change_minus_right(w, e_bias, average)) * (param.tunneling_exponent_right +
                                                                                                                                            param.tunneling_exponent_right * np.exp(energy_change_minus_right(w, e_bias, average) / param.kT) +
                                                                                                                                            (1 / param.kT)) * (fermi_distribution(energy_change_minus_right(w, e_bias, average)) ** 2)


def d_gamma_minus_right(w, e_bias, average=True):
    return param.gamma_zero_right * np.exp(param.tunneling_exponent_right * energy_change_minus_right(w, e_bias, average)) * (param.tunneling_exponent_right +
                                                                                                                     (param.tunneling_exponent_right - (1 / param.kT)) *
                                                                                                                     np.exp(energy_change_minus_right(w, e_bias, average) / param.kT)) * (fermi_distribution(energy_change_minus_right(w, e_bias, average)) ** 2)


# def d_gamma_plus_left(w, e_bias):
#     return 2 * param.gamma_zero_left * np.exp((param.tunneling_exponent_left + (1 / param.kT)) * energy_change_minus_left(w, e_bias)) * (param.tunneling_exponent_left +
#                                                                                                                                          param.tunneling_exponent_left * np.exp(energy_change_minus_left(w, e_bias) / param.kT) +
#                                                                                                                                          (1 / param.kT)) * (fermi_distribution(energy_change_minus_left(w, e_bias)) ** 2)
# def d_gamma_minus_left(w, e_bias):
#     return param.gamma_zero_left * np.exp(param.tunneling_exponent_left * energy_change_minus_left(w, e_bias)) * (param.tunneling_exponent_left +
#                                                                                                                   (param.tunneling_exponent_left - (1 / param.kT)) *
#                                                                                                                   np.exp(energy_change_minus_left(w, e_bias) / param.kT)) * (fermi_distribution(energy_change_minus_left(w, e_bias)) ** 2)


def gamma_plus(w, e_bias, epsilon_zero, width):
    return gamma_plus_left_modified(w, e_bias, epsilon_zero, width) + gamma_plus_right(w, e_bias)


def gamma_minus(w, e_bias, epsilon_zero, width):
    return gamma_minus_left_modified(w, e_bias, epsilon_zero, width) + gamma_minus_right(w, e_bias)


def d_gamma_plus(w, e_bias, epsilon_zero, width):
    return d_gamma_plus_left_modified(w, e_bias, epsilon_zero, width) + d_gamma_plus_right(w, e_bias)


def d_gamma_minus(w, e_bias, epsilon_zero, width):
    return d_gamma_minus_left_modified(w, e_bias, epsilon_zero, width) + d_gamma_minus_right(w, e_bias)


def gamma_total(w, e_bias, epsilon_zero, width):
    return gamma_plus(w, e_bias, epsilon_zero, width) + gamma_minus(w, e_bias, epsilon_zero, width)


def d_n(w, e_bias, epsilon_zero, width):
    return (1 / (gamma_total(w, e_bias, epsilon_zero, width) ** 2)) * (gamma_minus(w, e_bias, epsilon_zero, width) * d_gamma_plus(w, e_bias, epsilon_zero, width) - gamma_plus(w, e_bias, epsilon_zero, width) * d_gamma_minus(w, e_bias, epsilon_zero, width))


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


def kappa_tilde_mesh(w, e_bias, epsilon_zero, width):
    kappa_tilde_m = np.ndarray((len(w_list), len(e_bias_list)))

    for i in range(len(w_list)):
        for j in range(len(e_bias_list)):
            dn = d_n(w[i][j], e_bias[i][j], epsilon_zero, width)
            gamma_t = gamma_total(w[i][j], e_bias[i][j], epsilon_zero, width)

            _term = dn / gamma_t
            _k = ((np.cos(theta_list) ** 2) / np.pi) * (dn / gamma_t)

            k = np.trapz(_k, dx=d_theta)
            k *= (param.coupling_force * param.coupling_force / param.oscillator_mass)

            kappa_tilde_m[i][j] = k
    return kappa_tilde_m


def kappa_tilde_bias_slice(w, e_bias, epsilon_zero, width):
    kappa_tilde_b = np.ndarray(len(w_list))

    for i in range(len(w_list)):
        dn = d_n(w[i], e_bias, epsilon_zero, width)
        gamma_t = gamma_total(w[i], e_bias, epsilon_zero, width)

        _term = dn / gamma_t
        _k = ((np.cos(theta_list) ** 2) / np.pi) * (dn / gamma_t)

        k = np.trapz(_k, dx=d_theta)
        k *= (param.coupling_force * param.coupling_force / param.oscillator_mass)

        kappa_tilde_b[i] = k
    return kappa_tilde_b


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

    label_min = "(" f'{float(f"{minimum[0]:.3g}"):g}' + ", " + f'{float(f"{minimum[1]:.3g}"):g}' + ")"
    ax.plot(minimum[0], minimum[1], "wp")
    # ax.annotate(label, (minimum[0], minimum[1]), textcoords="offset points", xytext=(1, -12), ha="center")

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(0, 10)

    ax.axhline(8, linestyle='solid', color='g', linewidth=0.3)
    # points = ax.plot([-2, -1.3, -0.5, 1], [4, 4, 4, 4], 'gs')
    x, y = -2.4, 8
    point_label = "(" + str(x) + ", " + str(y) + ")"
    point = ax.plot(x, y, 'gx')

    handles = [mlines.Line2D([], [], color='w', marker='p', markersize=6, linewidth=0, markeredgecolor='black'), mlines.Line2D([], [], color='g', marker='x', markersize=6, linewidth=0)]
    ax.legend(handles, [label_min, point_label], loc='lower right')

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

    filename = "KappaDensity " + f"E={mechanical_energy}Wc. Min. Marked " + point_label + ".png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    fig.show()


def plot_kappa_tilde_symmetric():
    w_list2 = np.linspace(-5, 5, 200)
    e_bias_list2 = np.linspace(-10, 10, 200)
    w_mesh2, bias_mesh2 = np.meshgrid(w_list2, e_bias_list2)

    epsilon_zero = 0
    width = 0.5

    kt2 = kappa_tilde_mesh(w_mesh2, bias_mesh2, epsilon_zero, width)

    fig, ax = plt.subplots(num="Modified Left Rates Kappa")

    # ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='RdGy', alpha=1)
    # im = ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='coolwarm', alpha=1)
    # ax.axis(aspect='image')

    divnorm = colors.TwoSlopeNorm(vcenter=0)

    im = ax.contourf(w_mesh2, bias_mesh2, kt2, 20, cmap='plasma', norm=divnorm)  # Maybe plasma, viridis
    ax.contour(w_mesh2, bias_mesh2, kt2, 1, colors='black', levels=[0], linestyles='dashed')

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(-10, 10)

    ax.tick_params(top=True, right=True)
    ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
    ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

    title = "Modified Left Rates\n" + "E=" + str(mechanical_energy) + "$W_c$, $\epsilon_0=$" + str(epsilon_zero) + "$W_c$, $\sigma=$" + str(width) + "$W_c$, $A$=1"
    ax.set_title(title)

    cb = fig.colorbar(im, orientation='vertical')
    cb.set_label("$\kappa$-tilde / $s^{-1}$", size=14)
    # cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    ax.set_aspect(0.4)

    fig.tight_layout()
    kt = kappa_tilde_bias_slice(w_list, 6)
    # kt = kt / np.linalg.norm(kt)
    dn = d_n(w_list, 6)
    # dn = dn / np.linalg.norm(dn)

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

    fig.tight_layout()

# def plot_kappa_tilde_dn_dw_bias_slice():
#     kt = kappa_tilde_bias_slice(w_list, 6)
#     # kt = kt / np.linalg.norm(kt)
#     dn = d_n(w_list, 6)
#     # dn = dn / np.linalg.norm(dn)
#
#     fig, ax1 = plt.subplots()
#
#     color = 'tab:red'
#     ax1.set_xlabel('$W$ (units of $W_c$)', fontsize=14)
#     ax1.set_ylabel('$\kappa$-tilde ($s^{-1}$)', color=color, fontsize=14)
#     ax1.plot(w_list, kt, color=color)
#     ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
#
#     ax2 = ax1.twinx()
#
#     color = 'tab:blue'
#     ax2.set_ylabel('$\partial_W n$ ($W_c^{-1}$)', color=color, fontsize=14)
#     ax2.plot(w_list, dn, color=color)
#     ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
#
#     ax1.set_title("Comparing $\partial_W n$ and $\kappa$-tilde for $e V_b = 6 W_c$", fontsize=14)
#     plt.tight_layout()
#     fig.show()


# def plot_kappa_tildes():
#     global mechanical_energy
#
#     x, y = -0.5, 2.2
#     for i in range(4):
#         mechanical_energy = (i + 0) * 20
#
#         kt = kappa_tilde_mesh(w_mesh, bias_mesh)
#
#         fig, ax = plt.subplots()
#
#         # ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='RdGy', alpha=1)
#         # im = ax.imshow(dn, extent=[-5, 5, 0, 10], origin='lower', cmap='coolwarm', alpha=1)
#         # ax.axis(aspect='image')
#
#         im = ax.contourf(w_mesh, bias_mesh, kt, 20, cmap='coolwarm')
#         ax.contour(w_mesh, bias_mesh, kt, 1, colors='black', levels=[0], linestyles='dashed')
#
#         ax.plot(w_list, e_bias_left(w_list), "k")
#         ax.plot(w_list, e_bias_right(w_list), "k")
#         plt.ylim(0, 10)
#
#         point = ax.plot(x, y, 'mx')
#         handles = [mlines.Line2D([], [], color='m', marker='x', markersize=6, linewidth=0)]
#         ax.legend(handles, ["(" + str(x) + ", " + str(y) + ")"], loc='lower right')
#
#         ax.tick_params(top=True, right=True)
#         ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
#         ax.tick_params(axis='y', direction='in', length=6, labelsize=12)
#
#         ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
#         ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)
#
#         ax.set_title("Mechanical Energy: " + str(mechanical_energy) + "$W_c$")
#
#         cb = fig.colorbar(im, orientation='vertical')
#         cb.set_label("$\kappa$-tilde", size=14)
#         # cb.ax.set_title("$\partial_W n$", size=14)
#         cb.ax.tick_params(labelsize='large')
#
#         ax.set_aspect(0.75)
#         fig.show()


def plot_kappa_tilde_region_iii():
    global mechanical_energy

    epsilon_zero = 0
    width = 2

    mechanical_energy = 0
    kt0 = kappa_tilde_mesh(w_mesh, bias_mesh, epsilon_zero, width)

    _kt = np.sign(kt0)

    sign_changes = np.zeros((len(w_list), len(e_bias_list)))

    count = 3
    d_energy = 20
    for i in range(count):
        mechanical_energy = (i + 1) * d_energy

        kt = np.sign(kappa_tilde_mesh(w_mesh, bias_mesh, epsilon_zero, width))
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

    ax.axhline(6, linestyle='solid', color='g', linewidth=0.3)
    point = ax.plot([-2.7, -2.3, -1.57, 1.6], [6, 6, 6, 6], 'gs')

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

    E = count * d_energy
    filename = "RegioniiiKappaDensity " + f"Emax={E}Wc.png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    fig.show()


def plot_kappa_tilde_region_iii_symmetric():
    global mechanical_energy

    epsilon_zero = 0
    width = 2

    w_list2 = np.linspace(-5, 5, 200)
    e_bias_list2 = np.linspace(-10, 10, 200)
    w_mesh2, bias_mesh2 = np.meshgrid(w_list2, e_bias_list2)

    mechanical_energy = 0
    kt0 = kappa_tilde_mesh(w_mesh2, bias_mesh2, epsilon_zero, width)

    _kt = np.sign(kt0)

    sign_changes = np.zeros((len(w_list2), len(e_bias_list2)))

    count = 50
    d_energy = 50
    for i in range(count):
        mechanical_energy = (i + 1) * d_energy

        kt = np.sign(kappa_tilde_mesh(w_mesh2, bias_mesh2, epsilon_zero, width))
        for j in range(len(w_list2)):
            for k in range(len(e_bias_list2)):
                if _kt[j][k] != kt[j][k]:
                    sign_changes[j][k] += 1
        _kt = kt
    for j in range(len(w_list2)):
        for k in range(len(e_bias_list2)):
            if sign_changes[j][k] == 0:
                continue
            if sign_changes[j][k] == 1:
                sign_changes[j][k] = 0

    blue_colors = [(0, 0, 1, c) for c in np.linspace(0, 1, 100)]
    cmapblue = colors.LinearSegmentedColormap.from_list('mycmap', blue_colors, N=(np.amax(sign_changes) + 1))

    fig, ax = plt.subplots(num='Region iii')

    divnorm = colors.TwoSlopeNorm(vcenter=0)

    im = ax.contourf(w_mesh2, bias_mesh2, kt0, 20, cmap='plasma', norm=divnorm)  # Maybe plasma, viridis
    ax.contour(w_mesh2, bias_mesh2, kt0, 1, colors='black', levels=[0], linestyles='dashed')

    regions = ax.contourf(w_mesh2, bias_mesh2, sign_changes, 20, cmap=cmapblue)

    # ax.axhline(6, linestyle='solid', color='g', linewidth=0.3)
    # point = ax.plot([-2.7, -2.3, -1.57, 1.6], [6, 6, 6, 6], 'gs')

    ax.plot(w_list, e_bias_left(w_list), "k")
    ax.plot(w_list, e_bias_right(w_list), "k")
    plt.ylim(-10, 10)

    # handles = [mlines.Line2D([], [], color='m', marker='x', markersize=6, linewidth=0)]
    # ax.legend(handles, ["(" + str(x) + ", " + str(y) + ")"], loc='lower right')

    ax.tick_params(top=True, right=True)
    ax.tick_params(axis='x', direction='in', length=6, labelsize=12)
    ax.tick_params(axis='y', direction='in', length=6, labelsize=12)

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$e V_b$ (units of $W_c$)", fontsize=14)

    # ax.set_title("Mechanical Energy: " + str(mechanical_energy) + "$W_c$")

    # cb_regions = fig.colorbar(regions, ticks=np.linspace(0, 5, 6))
    # cb_regions.set_label("Roots")
    # cb_regions.ax.tick_params(labelsize='large')

    cb = fig.colorbar(im, orientation='vertical')
    cb.set_label("$\kappa$-tilde ($s^{-1}$)", size=14)
    # cb.ax.set_title("$\partial_W n$", size=14)
    cb.ax.tick_params(labelsize='large')

    title = "Max(E)= " + str(count * d_energy) + "$W_c$, $\epsilon_0=$" + str(epsilon_zero) + "$W_c$, $\sigma=$" + str(width) + "$W_c$"
    ax.set_title(title)

    ax.set_aspect(0.4)

    E = count * d_energy
    filename = "RegioniiiKappaDensity_" + f"Emax={E}Wc_symmetric_newgauss_left.png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    fig.tight_layout()


def plot_kappa_bias_slice_region_iii(e_bias):
    global mechanical_energy

    mechanical_energy = 0

    epsilon_zero = 0
    width = 2

    kt0 = kappa_tilde_bias_slice(w_list, e_bias, epsilon_zero, width)
    kt0_signs = np.sign(kt0)

    _kt_signs = np.sign(kt0)
    sign_changes = np.zeros(len(w_list))

    count = 50
    d_energy = 50

    start_iii_shade = None
    end_iii_shade = None
    start_ii_shade = None
    end_ii_shade = None
    for i in range(count):
        mechanical_energy = (i + 1) * d_energy

        kt = np.sign(kappa_tilde_bias_slice(w_list, e_bias, epsilon_zero, width))
        for j in range(len(w_list)):
            if _kt_signs[j] != kt[j]:
                sign_changes[j] += 1
        _kt_signs = kt
    for j in range(len(w_list)):
        if kt0_signs[j] == 0 and start_ii_shade is None:
            start_ii_shade = w_list[j]
        if kt0_signs[j] == 0 and start_ii_shade is not None:
            end_ii_shade = w_list[j]
        if kt0_signs[j] == -1 and kt0_signs[j] != kt0_signs[j - 1] and start_ii_shade is None:
            start_ii_shade = w_list[j]
        if kt0_signs[j] == 1 and kt0_signs[j] != kt0_signs[j - 1] and end_ii_shade is None:
            end_ii_shade = w_list[j]
        if sign_changes[j] == 0:
            continue
        if sign_changes[j] == 1:
            sign_changes[j] = 0
        if sign_changes[j] == 2 and sign_changes[j] != sign_changes[j - 1]:
            start_iii_shade = w_list[j]
        if sign_changes[j] == 0 and sign_changes[j] != sign_changes[j - 1]:
            end_iii_shade = w_list[j]

    fig, ax = plt.subplots(num='2D Region iii')

    ax.axvspan(start_iii_shade, end_iii_shade, color='b', alpha=0.4, lw=0, label="Region iii")
    ax.axvspan(start_ii_shade, end_ii_shade, color='m', alpha=0.4, lw=0, label="Region ii")
    ax.axhline(y=0, color='k', label='_nolegend_', alpha=0.5, linewidth=0.5)

    ax.set_xlabel('$W$ (units of $W_c$)', fontsize=14)
    ax.set_ylabel('$\kappa$-tilde ($s^{-1}$)', fontsize=14)
    ax.plot(w_list, kt0, color='k', label="$\kappa(E=0)$")


    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    ax.set_title("$E_{max}=$" + str(count * d_energy) + "$W_c$, $\sigma=$" + str(width) + "$W_c$, $\epsilon_0=$" + str(epsilon_zero) + "$W_c$, $eV_b=$" + str(e_bias) + "$W_c$", fontsize=14)

    ax.legend()

    fig.tight_layout()


def save_kappa_tilde_region_iii_data():
    global mechanical_energy

    mechanical_energy = 0
    kt0 = kappa_tilde_mesh(w_mesh, bias_mesh)

    _kt = np.sign(kt0)

    sign_changes = np.zeros((len(w_list), len(e_bias_list)))

    count = 60
    d_energy = 50
    for i in range(count):
        mechanical_energy = (i + 1) * d_energy

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

    ax.axhline(6, linestyle='solid', color='g', linewidth=0.3)
    point = ax.plot([-2.7, -2.3, -1.57, 1.6], [6, 6, 6, 6], 'gs')

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

    energy_max = count * d_energy
    data_name = "regioniii_data_" + f"Emax={energy_max}"
    np.save("../data/" + data_name, sign_changes)

    filename = "RegioniiiKappaDensity " + f"Emax={energy_max}Wc.png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    fig.show()


def plot_kappa_tilde_region_iii_from_data():
    global mechanical_energy

    mechanical_energy = 0
    kt0 = kappa_tilde_mesh(w_mesh, bias_mesh)

    sign_changes = np.load("../data/regioniii_data_Emax=3000.npy")

    blue_colors = [(0, 0, 1, c) for c in np.linspace(0, 1, 100)]
    cmapblue = colors.LinearSegmentedColormap.from_list('mycmap', blue_colors, N=(np.amax(sign_changes) + 1))

    fig, ax = plt.subplots(num='Region iii')

    divnorm = colors.TwoSlopeNorm(vcenter=0)

    im = ax.contourf(w_mesh, bias_mesh, kt0, 20, cmap='plasma', norm=divnorm)  # Maybe plasma, viridis
    ax.contour(w_mesh, bias_mesh, kt0, 1, colors='black', levels=[0], linestyles='dashed')

    regions = ax.contourf(w_mesh, bias_mesh, sign_changes, 20, cmap=cmapblue)

    ax.axhline(6, linestyle='solid', color='g', linewidth=0.3)

    xs = [-2.7, -2.3, -1.95, -1.57, 1.6, 2.8]
    ys = [6, 6, 6, 6, 6, 6]
    ns = ['a', 'b', 'c', 'd', 'e', 'f']
    points = ax.plot(xs, ys, 'gs')

    for i, txt in enumerate(ns):
        ax.annotate(txt, (xs[i], ys[i]), textcoords="offset points", xytext=(0, 6), ha="center", bbox=dict(boxstyle="circle,pad=0.05", fc=(1, 1, 1, 0.5), ec="none"))

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

    filename = "RegioniiiKappaDensity FromData3000.png"
    plt.savefig("../out/" + filename, bbox_inches='tight')
    print("Saved: " + filename)

    fig.show()


def plot_plus_rates_bias_slice(e_bias):
    epsilon_zero = 0
    width = 2

    gamma_left = np.ndarray(len(w_list))
    gamma_right = np.ndarray(len(w_list))

    for i in range(len(w_list)):
        gamma_left[i] = gamma_plus_left_modified(w_list[i], e_bias, epsilon_zero, width, False)
        gamma_right[i] = gamma_plus_right(w_list[i], e_bias, False)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), num='Modified Plus Rates Slice Subplots')

    axes[0].plot(w_list, gamma_left, 'g', label='$\Gamma^+_L$', linestyle='--')
    axes[1].plot(w_list, gamma_right, 'r', label='$\Gamma^+_R$', linestyle='-.')

    axes[0].xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    axes[1].xaxis.set_label_text("W (units of $W_c$)", fontsize=14)

    axes[0].yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    axes[0].tick_params(axis='x', direction='in', labelsize=12)
    axes[0].tick_params(axis='y', direction='in', labelsize=12)

    axes[1].tick_params(axis='x', direction='in', labelsize=12)
    axes[1].tick_params(axis='y', direction='in', labelsize=12)

    axes[0].legend()
    axes[1].legend()

    params = "$e V_b$=" + str(e_bias) + "$W_c$, $\epsilon_0$=" + str(epsilon_zero) + "$W_c$, $\sigma$=" + str(width) + "$W_c$"
    title = "Modified Tunneling Rates onto the Island\n" + params
    fig.suptitle(title)

    fig2, ax = plt.subplots(num='Modified Plus Rates Slice Single')

    ax.plot(w_list, gamma_left, 'g', label='$\Gamma^+_L$', linestyle='--')
    ax.plot(w_list, gamma_right, 'r', label='$\Gamma^+_R$', linestyle='-.')

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    ax.legend()

    ax.set_title(title)

    fig.tight_layout()
    fig2.tight_layout()


def plot_minus_rates_bias_slice(e_bias):
    epsilon_zero = 0
    width = 2

    gamma_left = np.ndarray(len(w_list))
    gamma_right = np.ndarray(len(w_list))

    for i in range(len(w_list)):
        gamma_left[i] = gamma_minus_left_modified(w_list[i], e_bias, epsilon_zero, width, False)
        gamma_right[i] = gamma_minus_right(w_list[i], e_bias, False)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), num='Modified Minus Rates Slice Subplots')

    axes[0].plot(w_list, gamma_left, 'g', label='$\Gamma^-_L$', linestyle='--')
    axes[1].plot(w_list, gamma_right, 'r', label='$\Gamma^-_R$', linestyle='-.')

    axes[0].xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    axes[1].xaxis.set_label_text("W (units of $W_c$)", fontsize=14)

    axes[0].yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    axes[0].tick_params(axis='x', direction='in', labelsize=12)
    axes[0].tick_params(axis='y', direction='in', labelsize=12)

    axes[1].tick_params(axis='x', direction='in', labelsize=12)
    axes[1].tick_params(axis='y', direction='in', labelsize=12)

    axes[0].legend()
    axes[1].legend()

    params = "$e V_b$=" + str(e_bias) + "$W_c$, $\epsilon_0$=" + str(epsilon_zero) + "$W_c$, $\sigma$=" + str(width) + "$W_c$"
    title = "Modified Tunneling Rates off of the Island\n" + params
    fig.suptitle(title)

    fig2, ax = plt.subplots(num='Modified Minus Rates Slice Single')

    ax.plot(w_list, gamma_left, 'g', label='$\Gamma^-_L$', linestyle='--')
    ax.plot(w_list, gamma_right, 'r', label='$\Gamma^-_R$', linestyle='-.')

    ax.xaxis.set_label_text("W (units of $W_c$)", fontsize=14)
    ax.yaxis.set_label_text("$\Gamma$ / Hz", fontsize=14)

    ax.legend()

    ax.set_title(title)

    fig.tight_layout()
    fig2.tight_layout()


# plot_plus_rates_bias_slice(5)
plot_kappa_bias_slice_region_iii(-4)
