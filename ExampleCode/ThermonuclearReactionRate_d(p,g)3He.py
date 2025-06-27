#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
# constants and definitions
M0 = 1.007276466583 
M1 = 2.01355319820		# masses (amu) of p and d
Z0 = 1
Z1 = 1				# charges of p and d

k_B = 0.086173324 # MeV/ GK #1.380649e-23 #J/K     
N_A = 6.02214076e23       

muc2 = 931.494061e0         # m_u*c^2 in MeV
clight = 299792458e2        # speed of light in cm/s

# reduced mass
mue = (M0 * M1)/(M0 + M1)

def Spoly(E,a,b,c,d):
    """
    Polynomial fit for the cross section
    """
    return a + E*b + (E**2)*c + (E**3)*d

def dpieta(E):
    return 0.98951013 * Z0 * Z1 * np.sqrt(mue/E)


# S-factor in MeV*b
S0 = 2.19e-7
S1 = 5.80e-6
S2 = 6.34e-6
S3 = -2.20e-6

                            
def integrand(E, T9):
    """
    integrand for the integral
    """
    return np.exp(-dpieta(E)) * Spoly(E, S0, S1, S2, S3) * np.exp(-E/(k_B*T9))

#def Nasv(T):
    """
 #   Reaction rate
    """
  #  sqrt_term = np.sqrt(8 / (np.pi * mue * muc2))
   # avog = N_A
    #prefactor = sqrt_term * clight * avog * 1e-24 *  ((k_B * T)**(-3/2))
    #integral, _ = quad(integrand, 0, 2.0, args=(T,))
    #return prefactor * integral
def Nasv(T):
    T9 = T  # Explicitly treating T as T9 in GK
    sqrt_term = np.sqrt(8 / (np.pi * mue * muc2))
    avog = N_A
    prefactor = sqrt_term * clight * avog * 1e-24 * ((k_B * T9)**(-3/2))

    # Constants for E0 and ΔE
    A = mue  # Effective reduced mass in amu

    E0 = 0.1220 * ((Z0**2 * Z1**2 * A)**(1/3)) * (T9)**(2/3)
    delta_E = 0.2368 * ((Z0**2 * Z1**2 * A)**(1/6)) * (T9)**(5/6)

    E_min = max(0, E0 - 2 * delta_E)  # prevent negative energy
    E_max = E0 + 2 * delta_E

    integral, _ = quad(integrand, E_min, E_max, args=(T9,))
    return prefactor * integral

T_range_GK = np.linspace(1e-3, .03, 100)
rates = [Nasv(T) for T in T_range_GK]

# Plotting!!!
plt.figure(figsize=(15, 7))
plt.plot(T_range_GK, rates, color='darkred')
plt.xlabel('Temperature (K)')
plt.ylabel(r'N$_A\langle\sigma v\rangle$ (cm$^3$ mol$^{-1}$ s$^{-1}$)')
plt.title('Reaction Rate vs Temperature for D(p,γ)³He')
#plt.yscale('log')
plt.grid(True)
plt.tight_layout()
plt.show()


# In[ ]:




