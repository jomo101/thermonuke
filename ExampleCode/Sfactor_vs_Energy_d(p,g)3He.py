#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt

# S-factor (MeV*b)
S0 = 2.19e-7
S1 = 5.80e-6
S2 = 6.34e-6
S3 = -2.20e-6

def Spoly(E, a, b, c, d):
    return a + E*b + E**2 * c + E**3 * d

# Energy range in MeV (log spaced for clarity)
E_MeV = np.logspace(-4, 0, 500)  # 0.0001 MeV to 1 MeV

# Calculate S-factor 
S_values = Spoly(E_MeV, S0, S1, S2, S3)

# Plottting
plt.figure(figsize=(12, 6))
plt.plot(E_MeV, S_values, color='navy', label='S-factor Polynomial Fit')
plt.xscale('log')
plt.xlabel('Energy (MeV)')
plt.ylabel('S-factor (MeV*b)')
plt.title(r'S-factor vs Energy for D(p,$\gamma$)$^3$He')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()


# In[ ]:




