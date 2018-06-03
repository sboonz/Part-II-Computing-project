from IsingHeader import energy_stdev
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

min_N = 10
max_N = 60
no_N = 11
N_for_C = np.linspace(min_N, max_N, no_N).astype(int)    # the value of N used for measuring heat capacity
MinT = 1.5
MaxT = 3.5
TNo = 41
T = np.linspace(MinT, MaxT, TNo)
spacing = int((MaxT - MinT) / (TNo - 1))

Tc_array = np.zeros(no_N)
Tc_error = np.linspace(spacing, spacing, no_N)

for i in range(no_N):

    C = np.array([(energy_stdev(T[x], N_for_C[i])**2)/(T[x]**2) for x in range(len(T))])

    for j in range(len(C)):
        if C[j] == np.amax(C):
            Tc_array[i] = T[j]
            break
    label_N = "N = " + str(N_for_C[i])
    plt.plot(T, C, label=label_N)
    plt.legend(loc='best')
plt.axvline(x=Tc_array[no_N-1], linestyle='--')
plt.xlabel(r'$T / J k_B^{{ -1 }}$')
plt.ylabel(r'$C_V / k_B$')
plt.title(r'Heat capacity plot to measure $T_C$')
SaveTitleC = 'CPlot.pdf'
plt.savefig("Figures/" + SaveTitleC)
plt.show()


def finite_scale_plot():

    # now fit a curve to determine T_c(infinity)

    def fit_func(n_in, tinf, alpha_fit, v):
        return tinf + alpha_fit * n_in ** (-1 / v)


    fit_params, fit_covariance = curve_fit(fit_func, N_for_C, Tc_array)  # the parameters are T_c(infinity), a and nu

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

    plt.plot(N_for_C, Tc_array, '+')
    plt.errorbar(N_for_C, Tc_array, yerr=Tc_error, fmt='b+')
    fit_label = r'${:01.3f} + {:01.3f} + N^{{\frac{{-1}}{{ {:01.3f} }} }} $'.format(T_c_infinity, alpha, nu)
    plt.plot(N_fine, T_c_fine, 'k', label=fit_label)
    plt.legend(loc='best')
    plt.plot(N_fine, T_c_fine, 'k', label=fit_label)
    plt.axhline(y=T_c_infinity, color='r', linestyle='--')
    plt.xlabel("N")
    plt.ylabel(r"$T_C / J k_B^{{ -1 }}$")
    plt.title("Finite-size scaling plot for critical temperature from Tc plot")
    plt.savefig("Figures/FiniteScaling/PlotCV")
    plt.show()

    f = open("ParametersFiniteScalingTcCV.txt", "w")
    f.write("T_c(infinity) = " + str(T_c_infinity) + " +/- " + str(T_c_infinity_error) + "\n")
    f.write("a = " + str(alpha) + " +/- " + str(alpha_error) + "\n")
    f.write("nu = " + str(nu) + " +/- " + str(nu_error))
    f.close()

finite_scale_plot()