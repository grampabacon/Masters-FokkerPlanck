# Oscillator variables
oscillator_mass = 21.3e-21  # kg
oscillator_frequency = 10e9  # Hertz
quality_factor = 10e5

# System variables
coupling_parameter = 170
coupling_force = 19.54123e-12  # 6.179485810162524e-13
temperature = 8  # Kelvin

# Energy scale
W_c = 5e-3 * 1.60217662e-19

# Voltages
bias_voltage = (1 * W_c) / 1.60217662e-19
gate_voltage = 0

# Capacitance
left_capacitance = 0.3 * 10e-15
right_capacitance = 0.75 * 10e-15

# Energies
W = W_c
W_L = -1 * 1.60217662e-19 * (bias_voltage / 2)
W_R = 1 * 1.60217662e-19 * (bias_voltage / 2)

# Tunneling rates
gamma_zero_left = 125e9  # Hertz
gamma_zero_right = 125e9  # Hertz
tunneling_exponent_left = 0.3
tunneling_exponent_right = 0.75
