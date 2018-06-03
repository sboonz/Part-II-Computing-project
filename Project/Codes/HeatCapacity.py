from IsingHeader import energy_stdev
import numpy as np
import matplotlib.pyplot as plt


N_for_C = 16    # the value of N used for measuring heat capacity
MinT = 1
MaxT = 6
TNo = 100
T = np.linspace(MinT, MaxT, TNo)

CSamples = 5   # number of heat capacity values to average over
CMatrix = np.array([[(energy_stdev(T[y], N_for_C)**2)/(T[y]**2) for x in range(CSamples)] for y in range(len(T))])
C = np.array([np.mean(CMatrix[x]) for x in range(len(T))])
C_error = np.array([np.std(CMatrix[x]) for x in range(len(T))])

T_c_est = 0 # estimate for T_c_infinity

for i in range(len(C)):
    if C[i] == np.amax(C):
        T_c_est = T[i]
        break

plt.plot(T, C, 'k')
plt.axvline(x=T_c_est, color='r', linestyle='--')
plt.errorbar(T, C, yerr=C_error, fmt='b+')
plt.xlabel(r'$T / J k_B^{{ -1 }}$')
plt.ylabel(r'$C_V / k_B$')
plt.title(r'Heat capacity plot to measure $T_C$')
SaveTitleC = 'CPlot.pdf'
plt.savefig("Figures/" + SaveTitleC)


error_Tc = (MaxT - MinT) / TNo
f = open("TcFromHeatCapacity.txt", "w")
f.write("T_c = " + str(T_c_est) + " +/- " + str(error_Tc))
f.close()

