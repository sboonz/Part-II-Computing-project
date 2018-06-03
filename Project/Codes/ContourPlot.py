from IsingHeader import spin_plot
from IsingHeader import energy_plot
from IsingHeader import N_init

TArray = [0.1, 1, 2, 2.2, 2.5, 20]	# selected temperatures in J/K_B

for i in range(len(TArray)):
	spin_plot(TArray[i], N_init)
	energy_plot(TArray[i], N_init)
