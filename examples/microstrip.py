import numpy as np
import scipy.constants
from philsol.classy import phil_class
import matplotlib.pyplot as mplt


# create permittivity matrix of a microstrip cross section
def n_crosssection(x, y):
    # center of line
    x0 = 5
    y0 = 1

    # width and height
    w = 1
    h = 1.5

    # conductor thickness (choose resolution accordingly)
    d_cond = 0.05

    # substrate loss tangent and real relative permittivity
    tand_subst = 0.02
    eps_subst = 4.5 * (1 + 1j * tand_subst)

    # conductor conductivity
    eps_cond = 1e9j

    eps = np.ones((len(x), len(y), 3), dtype=np.complex)
    for j, y_j in enumerate(y):
        for i, x_i in enumerate(x):
            if -0.5 * h < y_j - y0 < 0.5 * h:
                eps[i, j] = [eps_subst, eps_subst, eps_subst]
            elif -0.5 * h - d_cond <= y_j - y0 <= -0.5 * h:
                if -20 < x_i - x0 < 20:
                    eps[i, j] = [eps_cond, eps_cond, eps_cond]
            elif 0.5 * h <= y_j - y0 <= 0.5 * h + d_cond:
                if -0.5 * w < x_i - x0 < 0.5 * w:
                    eps[i, j] = [eps_cond, eps_cond, eps_cond]
            else:
                eps[i, j] = [1.0, 1.0, 1.0]
    return np.sqrt(eps)


# create cross section
x = np.linspace(0, 10, 200)
y = np.linspace(0, 10, 200)

# calculate resolution from x and y arrays
dx = x[1] - x[0]
dy = y[1] - y[0]

# get matrix of refractive indices of cross section
n = n_crosssection(x, y)

# set frequency and layout unit
freq_Hz = 10e9
layout_unit_Meter = 1e-3

# calculate initial guess of propagation constant (in vacuum)
k0 = 2 * np.pi * freq_Hz / scipy.constants.c * layout_unit_Meter

# plot refractive profile
mplt.figure()
mplt.pcolormesh(np.abs(n[:, :, 0].transpose()))
mplt.show()

# create and start simulation
line = phil_class(n=n, k0=k0, dx=dx, dy=dy)
line.build_stuff(matrices=True)
line.solve_stuff(neigs=1, beta_trial=k0 * np.sqrt(4.5), extra_fields=True)

# reshape simulated E and H
E = np.reshape(line.E, (line.Eigs, line.num_x, line.num_y, 3))
H = np.reshape(line.H, (line.Eigs, line.num_x, line.num_y, 3))

# convert Gaussian to SI units (correct?)
E = E / np.sqrt(4 * np.pi * scipy.constants.epsilon_0)
H = H / np.sqrt(4 * np.pi * scipy.constants.mu_0)

# plot Ex, Ey, Ez in cross section
mplt.figure()
mplt.pcolormesh(np.abs(E[0, :, :, 0]))
mplt.show()
mplt.pcolormesh(np.abs(E[0, :, :, 1]))
mplt.show()
mplt.pcolormesh(np.abs(E[0, :, :, 2]))
mplt.show()

# calculate and print line properties
k = line.beta[0]
att = np.exp(np.imag(line.beta[0]))
att_per_meter = att / layout_unit_Meter
att_per_meter_db = 20 * np.log10(att) / layout_unit_Meter
wavelength_meter = 2 * np.pi / np.real(k) * layout_unit_Meter
eps_r_eff = (scipy.constants.c / freq_Hz / wavelength_meter) ** 2

print('wave vector: k = {}'.format(k))
print('attenuation: {} dB/m'.format(att_per_meter_db))
print('guided wavelength: {} m'.format(wavelength_meter))
print('effective relative permittivity: {}'.format(eps_r_eff))
