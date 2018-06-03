from IsingHeader import auto_corr
from IsingHeader import tau_array_gen
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import spline

''' for this task, we will manually change the step size to about 50, since for large steps we can't see the decay easily
and by reducing the steps, we reduce the run time '''
min_T = 1.25	# temperatures are near critical
max_T = 4.25
no_T = 4
TArray = np.linspace(min_T, max_T, no_T)	# temperatures for the plot
min_N = 15
max_N = 60
no_N = 4
N_values_float = np.linspace(min_N, max_N, no_N)
N_values = N_values_float.astype(int)

# this will be needed for the tau_e against T plot:

no_T_fine = 40
TArray_fine = np.linspace(min_T, max_T, no_T_fine)	# temperatures for the plot

# tau array

tau = tau_array_gen()	# tau array is imported from the 'header' file
taue_line = np.linspace(np.exp(-1), np.exp(-1), len(tau))

# initializing a matrix for tau_e values for different T and N

taue_tuple = (no_T_fine, no_N)
taue_array = np.zeros(taue_tuple)

# function to find the point where tau_e is

def find_taue(a_in):
	taue_value = 0
	for i in range(len(a_in)):
		if a_in[i] <= np.exp(-1):
			taue_value = tau[i]
			break
	return taue_value


for i in range(no_N):
	plt.subplot(no_N, 1, i+1)
	plot_title = "N = " + str(N_values[i])
	plt.title(plot_title)
	for j in range(no_T):
		a_mean = auto_corr(TArray[j], N_values[i])
		plt.plot(tau, a_mean, label=str(TArray[j]))
		plt.plot(tau, taue_line,"--")
		plt.legend(loc='best')
plt.xlabel(r'No of delay steps, $ \tau $, in arbitrary units')
plt.ylabel('a')
plt.savefig("Figures/AutoCorrelationOverTime.pdf")
plt.show()

taue_array = np.array([[find_taue(auto_corr(TArray_fine[x], N_values[y])) for x in range(no_T_fine)] for y in range(no_N)])

for i in range(no_N):
	y = taue_array[i,:]	# created a dummy variable for smoothing
	plt.subplot(no_N, 1, i+1)
	plt.plot(TArray_fine, y, label=str(N_values[i]))
	plt.legend(loc='best')
plt.xlabel(r"$T_C / J k_B^{{ -1 }}$")
plt.ylabel(r'$ \tau_e $ in steps')
plot_title = r'Plot of $ \tau_e $ against T for different N'
plt.title(plot_title)
plt.savefig("Figures/TaueAgainstTPlot.pdf")
plt.show()

