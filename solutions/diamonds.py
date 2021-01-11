import main
import numpy as np
import inputs.constants as const
import inputs.parameters as param
import matplotlib.pyplot as plt

#############################################################
#                                                           #
# See stationary_solution.py for the Fokker-Planck solution #
#                                                           #
#############################################################

# https://www.tutorialspoint.com/matplotlib/matplotlib_contour_plot.htm

# # Bias voltage range
# bias_min = - 5e-3 * param.W_c / const.ELEMENTARY_CHARGE
# bias_max = 10e-3 * param.W_c / const.ELEMENTARY_CHARGE
# bias_divisions = 1000
# bias_list = np.linspace(bias_min, bias_max, bias_divisions)
# d_bias = (bias_max - bias_min) / bias_divisions

# W range
w_min = - 5 * param.W_c
w_max = 5 * param.W_c
w_divisions = 1000
w_list = np.linspace(w_min, w_max, w_divisions)
d_w = (w_max - w_min) / w_divisions

gate_capacitance = 0
left_capacitance = 10e-15
right_capacitance = 10e-15


def total_capacitance():
    return gate_capacitance + left_capacitance + right_capacitance


def e_v_b_left():
    return (total_capacitance() / right_capacitance) * ((const.ELEMENTARY_CHARGE ** 2 / (2 * total_capacitance())) - w_list)


def e_v_b_right():
    return (total_capacitance() / left_capacitance) * (- (const.ELEMENTARY_CHARGE ** 2 / (2 * total_capacitance())) + w_list)


plt.plot(w_list / param.W_c, e_v_b_left() / param.W_c, "k")
plt.plot(w_list / param.W_c, e_v_b_right() / param.W_c, "g")
plt.xlabel("W / W_c")
plt.ylabel("e V_b / W_c")
plt.show()
