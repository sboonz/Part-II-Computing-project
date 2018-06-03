import numpy as np
import matplotlib.pyplot as plt
from IsingHeader import magnetisation_at_equilibrium
from IsingHeader import energy_at_equilibrium
from scipy.optimize import curve_fit

# temperature array

MinT = 2.2
MaxT = 2.6
no_of_data = 21
TArray = np.linspace(MinT, MaxT, no_of_data)
spacing = (MaxT - MinT) / (no_of_data - 1)

# N array

min_N = 15
max_N = 55
No_N = 9
N_values_float = np.linspace(min_N, max_N, No_N)
N_values = N_values_float.astype(int)

# matrices for the plots: T for rows, N for columns
magnetisation_matrix_eq = np.array([[magnetisation_at_equilibrium(TArray[x], N_values[y]) for x in range(no_of_data)] for y in range(No_N)])

energy_matrix_eq = np.array([[energy_at_equilibrium(TArray[x], N_values[y]) for x in range(no_of_data)] for y in range(No_N)])

# finding the value of Tc

half_mark = 0.5

''' since magnetisation is normalized, we can take the value for when equilibrium magnetisation is 1/2
to be the transition temperature- since free energy change = 0 => both perfectly disordered and ordered states are 
equally likely so averages to 1/2 '''

half_line = np.array([0.5 for x in range(no_of_data)])  # to determine transition point

Tc_array = np.zeros(No_N)

Tc_error_array = np.zeros(No_N)


def Fermi_Dirac_fit(t, alpha, beta):
    return 1/(np.exp((t-alpha)/beta) + 1)

def Tc_value_finder(input_array):
    FD_param, FD_covar = curve_fit(Fermi_Dirac_fit, TArray, input_array)
    return FD_param[0], FD_param[1], FD_covar[0][0]**0.5    # key: the outputs are alpha, beta and error in alpha i.e Tc

def plot_M_T():
    for i in range(No_N):
        mag_array_for_plot = magnetisation_matrix_eq[i, :]

        # plot parameters:
        parameters_list = Tc_value_finder(mag_array_for_plot)
        alpha = parameters_list[0]
        beta = parameters_list[1]
        Fermi_Dirac_fit = 1/(np.exp((TArray - alpha)/beta) + 1)

        # assign Tc values
        Tc_array[i] = alpha
        Tc_error_array[i] = parameters_list[2]

        # plot of the magnetisation
        plt.subplot(No_N, 1, i+1)
        plt.plot(TArray, mag_array_for_plot, 'k+', TArray, Fermi_Dirac_fit, 'b', TArray, half_line, 'k--')
        plt.axvline(x=alpha, color='r', linestyle='--')
        N_title = "N = " + str(N_values[i])
        plt.title(N_title)
        plt.xlabel(r"$T / J k_B^{{ -1 }}$")
        plt.ylabel("M")
    MTitle = 'Figures/MAainstTForN.pdf'
    plt.savefig(MTitle)
    plt.show()


def plot_E_T():
    for i in range(No_N):
        en_array_for_plot = energy_matrix_eq[i, :]
        plt.subplot(No_N, 1, i+1)
        plt.plot(TArray, en_array_for_plot, 'k--')
        plt.axvline(x=Tc_array[i], color='r', linestyle='--')
        N_title = "N = " + str(N_values[i])
        plt.title(N_title)
        plt.xlabel(r"$T / J k_B^{{ -1 }}$")
        plt.ylabel("Mean energy / J")
    ETitle = 'Figures/EAgainstT.pdf'
    plt.savefig(ETitle)
    plt.show()


def finite_scale_plot():

    # now fit a curve to determine T_c(infinity)

    def fit_func(n_in, tinf, alpha_fit, v):
        return tinf + alpha_fit * n_in ** (-1 / v)


    fit_params, fit_covariance = curve_fit(fit_func, N_values, Tc_array)  # the parameters are T_c(infinity), a and nu

    T_c_infinity = fit_params[0]
    T_c_infinity_error = fit_covariance[0][0] ** 0.5
    alpha = fit_params[1]
    alpha_error = fit_covariance[1][1] ** 0.5
    nu = fit_params[2]
    nu_error = fit_covariance[2][2] ** 0.5

    # for the best fit, we want a smooth curve, so need a finer array

    sample_number = 256
    N_fine = np.linspace(min_N, max_N, 256)

    T_c_fine = T_c_infinity + alpha * N_fine ** (-1 / nu)

    plt.plot(N_values, Tc_array, '+')
    plt.errorbar(N_values, Tc_array, yerr=Tc_error_array, fmt='k+')
    fit_label = r'${:01.3f} + {:01.3f} N^{{\frac{{-1}}{{ {:01.3f} }} }} $'.format(T_c_infinity, alpha, nu)
    plt.plot(N_fine, T_c_fine, 'k', label=fit_label)
    plt.legend(loc='best')
    plt.axhline(y=T_c_infinity, color='r', linestyle='--')
    plt.xlabel("N")
    plt.ylabel(r"$T_C / J k_B^{{ -1 }}$")
    plt.title("Finite-size scaling plot for critical temperature from Tc plot")
    plt.savefig("Figures/FiniteScaling/Plot")
    plt.show()

    f = open("ParametersFiniteScalingTc.txt", "w")
    f.write("T_c(infinity) = " + str(T_c_infinity) + " +/- " + str(T_c_infinity_error) + "\n")
    f.write("a = " + str(alpha) + " +/- " + str(alpha_error) + "\n")
    f.write("nu = " + str(nu) + " +/- " + str(nu_error))
    f.close()

plot_M_T()
plot_E_T()
finite_scale_plot()