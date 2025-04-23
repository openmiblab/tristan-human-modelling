import numpy as np
import matplotlib.pyplot as plt

# Constants
S0 = 100
actualFA = 12 # deg
estimatedFA = 15 #deg
TR = 0.0037 #sec
R10 = 1 #Hz
Cmax = 5 #mM
r1 = 5 #Hz/mM

# Concentration and R1
C = np.arange(0, Cmax, 0.1)
R1 = R10 + r1*C

# Generate a signal with the actual FA
cFA = np.cos(actualFA*np.pi/180)
E = np.exp(-TR*R1)
Sactual = S0 * ( 1 - E) / (1 - cFA * E )

# Reconstruct the concentrations with an estimated FA
cFA = np.cos(estimatedFA*np.pi/180)
E0 = np.exp(-TR*R10)
F0 = ( 1 - E0) / (1 - cFA * E0 ) # !! Flip angle used here !!
Y = F0 * Sactual / Sactual[0]

# Y = (1 - E) / (1 - cFA * E )
# Solve for E:
Erecon = (1 - Y) / (1 - cFA * Y)  # !! Flip angle used here !!

# Derive conc:
R1recon = -np.log(Erecon)/TR
Crecon = (R1recon - R1recon[0]) / r1

# Show results
plt.plot(C, C, 'k--')
plt.plot(C, Crecon, 'r-')
plt.xlabel('Actual concentration (mM)')
plt.ylabel('Reconstructed concentration (mM)')
plt.show()




