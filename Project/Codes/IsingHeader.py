import random as rd
import numpy as np
import matplotlib.pyplot as plt

# physical parameters for magnetic moment energy
N_init = 50  # this is the default value for N
J = 1  # in units of KB K i.e. J = 1 = 1.38 x 10^23 J -> we have set KB = 1 for our 'natural units'
mu = 1  # assuming non-magnetic medium
H = 0

# time array

time_steps = 200  # from the magnetisation graphs for low temperatures, we see equilibrium is
# established after around 20-30 steps, so need to have max no. of steps bigger than that
TimeArray = np.linspace(0, time_steps, time_steps)

equilibrium_step = 100  # rough estimate for no. of steps to eq. from M plot for MOST values of T


# generate a matrix representing the spins

def s_gen(N):
    return np.array([[1 for x in range(N)] for y in range(N)]) # s for spin lattice


# function to give the energy for a particular spin matrix

def e_adj(X, i, j, N):
    return -J * X[i][j] * ((X[i][(j + 1) %N] + X[i][(j - 1)%N] + X[(i + 1)%N][j] + X[(i - 1)%N][j]) + mu * H)


def delta_e(X, i, j, N):
    return -2 * e_adj(X, i, j, N)


# if the spin is opposite of what it is, then the energy should be opposite to the origin orientation
# hence delta_E = 2E, where E is defined by the equation for spin energy in the Ising model

# create an array to see the change in magnetisation overtime

# spin matrix change after a time step

def change_spin_function(X, Tin, N):  # X is a particular spin matrix, Tin is the temperature
    Xc = X  # updates the matrix every time a spin is flipped, which will affect
    # the energy of the neighbours
    for i in range(N):
        for j in range(N):
            p = rd.uniform(0, 1)
            de = delta_e(Xc, i, j, N)
            if de < 0 or np.exp(-de / Tin) > p:
                Xc[i][j] = -Xc[i][j]
            else:
                Xc[i][j] = Xc[i][j]
    return Xc

# magnetisation- net spin of the lattice


def magnetisation(X):
    return abs(np.mean(X))


def mean_energy(X):  # here X needs to be a 2D energy matrix
    return np.mean(X)


def matrix_generator(Tin, N):
    """
    Assign values to spin and energy matrices for a given temp. Tin
    """
    array_tuples = (time_steps, N, N)
    spin_dummy = np.zeros(array_tuples)
    energy_dummy = np.zeros(array_tuples)
    # spin matrix:
    spin_dummy[0] = s_gen(N)
    for i in range(1, time_steps):
        spin_dummy[i, :, :] = change_spin_function(spin_dummy[i - 1, :, :], Tin, N)

    # energy matrix:
    for i in range(time_steps):
        energy_dummy[i, :, :] = np.array([[e_adj(spin_dummy[i], x, y, N) for x in range(N)] for y in range(N)])
    return (spin_dummy, energy_dummy)


# a giant function/method for plotting magnetisation against time and spin contour

def spin_plot(Tin, N):
    spin_over_time = matrix_generator(Tin, N)[0]
    mag_array = np.array([magnetisation(spin_over_time[x]) for x in range(time_steps)])
    plt.plot(TimeArray, mag_array)
    mag_title = 'Magnetisation over time for T = ' + str(Tin) + ' J/K_B'
    plt.title(mag_title)
    plt.xlabel('time steps / arb. units')
    plt.ylabel('Magnetisation / (J/ spin)')
    save_title = 'MagnetisationT' + str(Tin) + 'N' + str(N) + '.pdf'
    plt.savefig('Figures/' + save_title)
    plt.show()

    # plot of the spin at select times
    selected_steps = [1, 5, 20, 100]
    plot_number = len(selected_steps)
    contour_title = 'Evolution of spin lattice over time for N = ' + str(N) + ', T = ' +str(Tin)
    plt.title(contour_title)
    for i in range(plot_number):
        plt.subplot(2, 2, i + 1)
        plt.pcolormesh(spin_over_time[selected_steps[i] - 1], vmin=-1, vmax=1)
        title = 'After ' + str(selected_steps[i]) + ' step'
        plt.title(title)
        plt.colorbar()
    save_title2 = 'SpinOverTimeT' + str(Tin) + 'N' + str(N) + '.pdf'
    plt.savefig('Figures/' + save_title2)
    plt.show()


# generate mean equilibrium magnetisation values for M against T plot

def magnetisation_at_equilibrium(Tin, N):  # returns the mean magnetisation at equilibrium
    spin_over_time = matrix_generator(Tin, N)[0]
    mag_array = np.array([magnetisation(spin_over_time[x]) for x in range(equilibrium_step, time_steps)])
    return np.mean(mag_array)


# method for plotting energy and its contour against time

def energy_plot(Tin, N):
    energy_matrix = matrix_generator(Tin, N)[1]
    energy_array = np.array([mean_energy(energy_matrix[x]) for x in range(time_steps)])
    plt.plot(TimeArray, energy_array)
    energy_title = 'Mean energy over time for N = ' + str(N) + ', T = ' + str(Tin)
    plt.title(energy_title)
    plt.xlabel('No of time steps')
    plt.ylabel('Mean energy in unit of J')
    save_title = 'MeanEnergyT' + str(Tin) + 'N' + str(N) + '.pdf'
    plt.savefig('Figures/' + save_title)
    plt.show()

    # plot of the energy of the lattice at select times

    selected_steps = [1, 5, 20, 100]
    plot_number = len(selected_steps)
    contour_title = 'Evolution of energy distribution over time for N = ' + str(N) + ', T = ' + str(Tin)
    plt.title(contour_title)
    for i in range(plot_number):
        plt.subplot(2, 2, i + 1)
        plt.pcolormesh(energy_matrix[selected_steps[i] - 1], vmin=-4, vmax=4)
        title = 'After ' + str(selected_steps[i]) + ' step'
        plt.title(title)
        plt.colorbar()
    save_title2 = 'EnergyOverTimeT' + str(Tin) + 'N' + str(N) + '.pdf'
    plt.savefig('Figures/' + save_title2)
    plt.show()


def energy_at_equilibrium(Tin, N):
    energy_matrix = matrix_generator(Tin, N)[1]
    energy_array = np.array([mean_energy(energy_matrix[x]) for x in range(equilibrium_step, time_steps)])
    return np.mean(energy_array)  # gives mean total energy


def energy_stdev(Tin, N):  # need this to measure heat capacity via fluctuation dissipation
    energy_matrix = matrix_generator(Tin, N)[1]
    energy_array = np.array([mean_energy(energy_matrix[x]) for x in range(equilibrium_step, time_steps)])
    return np.std(energy_array)


def auto_corr(Tin, N):   # generates array for auto-correlation
    # first two steps same as SpinMatrixGen
    spin_over_time = matrix_generator(Tin, N)[0]
    mag_array = np.array([magnetisation(spin_over_time[x]) for x in range(time_steps)])
    tau = np.linspace(0, time_steps - equilibrium_step, time_steps - equilibrium_step)
    mean_M = np.mean(mag_array[equilibrium_step: time_steps])
    # we need an array for the magnetisation after equilibrium is established
    equilibrium_mag = np.array([(mag_array[x] - mean_M) for x in range(equilibrium_step, time_steps)])
    # numpy built-in auto-correlation (cross-correlation) function
    Afull = np.correlate(equilibrium_mag, equilibrium_mag, mode='full')
    A = Afull[int(len(Afull) / 2):len(Afull)]
    a = abs(A / A[0])
    return a


def tau_array_gen():
    return np.linspace(0, time_steps - equilibrium_step, time_steps - equilibrium_step)

