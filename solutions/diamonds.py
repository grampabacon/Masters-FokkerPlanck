import main
import numpy as np
import inputs.constants as const
import inputs.parameters as param
import matplotlib.pyplot as plt

# https://www.tutorialspoint.com/matplotlib/matplotlib_contour_plot.htm

# Bias voltage range
bias_min = - 100e-3  # * param.W_c
bias_max = 100e-3  # * param.W_c
bias_divisions = 1000
bias_list = np.linspace(bias_min, bias_max, bias_divisions)
d_bias = (bias_max - bias_min) / bias_divisions

# Gate voltage range
gate_min = - 5e-3  # * param.W_c
gate_max = 5e-3  # * param.W_c
gate_divisions = 1000
gate_list = np.linspace(bias_min, bias_max, bias_divisions)
d_gate = (bias_max - bias_min) / bias_divisions

gate_capacitance = 10e-15
left_capacitance = 0.3 * 10e-15
right_capacitance = 0.75 * 10e-15


def total_capacitance():
    return gate_capacitance + left_capacitance + right_capacitance


def energy_change_left_plus(n, bias_voltage, gate_voltage):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / 2) + (n - (gate_voltage * gate_capacitance) / const.ELEMENTARY_CHARGE) - (right_capacitance * bias_voltage) / const.ELEMENTARY_CHARGE)


def energy_change_left_minus(n, bias_voltage, gate_voltage):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / 2) - (n - (gate_voltage * gate_capacitance) / const.ELEMENTARY_CHARGE) + (right_capacitance * bias_voltage) / const.ELEMENTARY_CHARGE)


def energy_change_right_plus(n, bias_voltage, gate_voltage):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / 2) + (n - (gate_voltage * gate_capacitance) / const.ELEMENTARY_CHARGE) + (left_capacitance * bias_voltage) / const.ELEMENTARY_CHARGE)


def energy_change_right_minus(n, bias_voltage, gate_voltage):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance()) * ((1 / 2) - (n - (gate_voltage * gate_capacitance) / const.ELEMENTARY_CHARGE) - (left_capacitance * bias_voltage) / const.ELEMENTARY_CHARGE)


elp_vec = np.vectorize(energy_change_left_plus)
elm_vec = np.vectorize(energy_change_left_minus)

erp_vec = np.vectorize(energy_change_right_plus)
erm_vec = np.vectorize(energy_change_right_minus)


x, y = np.meshgrid(bias_list, gate_list)
print(x)
#z = elp_vec(0, x, y)
z = erp_vec(0, x, y)

# plt.plot(bias_list / param.W_c, energy_change_left_plus(0, bias_list, gate_list), "g")
plt.contourf(x, y, z)

cb = plt.colorbar()
cb.set_label("Energy Change")

plt.xlabel("Bias voltage / W_c")
plt.ylabel("Gate voltage / W_c")
plt.show()
