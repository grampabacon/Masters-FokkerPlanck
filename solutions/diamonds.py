import main
import numpy as np
import inputs.constants as const
import inputs.parameters as param

# https://www.tutorialspoint.com/matplotlib/matplotlib_contour_plot.htm

# Bias voltage range
bias_min = - 5 * param.W_c
bias_max = 5 * param.W_c
bias_divisions = 1000
bias_list = np.linspace(bias_min, bias_max, bias_divisions)
d_bias = (bias_max - bias_min) / bias_divisions

# Gate voltage range
gate_min = - 5 * param.W_c
gate_max = 5 * param.W_c
gate_divisions = 1000
gate_list = np.linspace(bias_min, bias_max, bias_divisions)
d_gate = (bias_max - bias_min) / bias_divisions

# Theta for averaging over phase.
theta_divisions = 500
thetas = np.linspace(0, 2 * np.pi, theta_divisions)
d_theta = (2 * np.pi) / theta_divisions

# Range of mechanical energies for plotting
energy_min = 0 * param.W_c
energy_max = 1000 * param.W_c
energy_divisions = 1000
mechanical_energies = np.linspace(energy_min, energy_max, energy_divisions)
d_energy = (energy_max - energy_min) / energy_divisions

def gate_capacitance(mechanical_energy):
    return 10e-15 * (1)


def total_capacitance(mechanical_energy):
    return gate_capacitance(mechanical_energy) + param.left_capacitance + param.right_capacitance


def energy_change_plus_left(mechanical_energy):
    return (const.ELEMENTARY_CHARGE ** 2 / total_capacitance(mechanical_energy)) * ((1 / 2) + (main.average_occupation(mechanical_energy) + (gate_capacitance(mechanical_energy) * param.gate_voltage / const.ELEMENTARY_CHARGE)) + (
                param.right_capacitance * param.bias_voltage / const.ELEMENTARY_CHARGE)) - param.coupling_force * main.displacement(mechanical_energy)
