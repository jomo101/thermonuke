#!/usr/bin/env python
# coding: utf-8

# In[3]:


import gvar as gv

# Define gvars
S0_g = gv.gvar(2.19e-7, 0.10e-7)
S1_g = gv.gvar(5.80e-6, 0.24e-6)
S2_g = gv.gvar(6.34e-6, 0.88e-6)
S3_g = gv.gvar(-2.20e-6, 0.52e-6)

# Sampling function
def sample_S_params(n_samples=100):
    rng = np.random.default_rng()
    S0_samples = rng.uniform(S0_g.mean - S0_g.sdev, S0_g.mean + S0_g.sdev, n_samples)
    S1_samples = rng.uniform(S1_g.mean - S1_g.sdev, S1_g.mean + S1_g.sdev, n_samples)
    S2_samples = rng.uniform(S2_g.mean - S2_g.sdev, S2_g.mean + S2_g.sdev, n_samples)
    S3_samples = rng.uniform(S3_g.mean - S3_g.sdev, S3_g.mean + S3_g.sdev, n_samples)
    return list(zip(S0_samples, S1_samples, S2_samples, S3_samples))

def Nasv_with_S(T, Sparams):
    S0, S1, S2, S3 = Sparams
    T9 = T  # Make sure T is in GK
    sqrt_term = np.sqrt(8 / (np.pi * mue * muc2))
    prefactor = sqrt_term * clight * N_A * 1e-24 * ((k_B * T9)**(-3/2))
    
    A = mue  # effective reduced mass

    E0 = 0.1220 * ((Z0**2 * Z1**2 * A)**(1/3)) * (T9)**(2/3)
    delta_E = 0.2368 * ((Z0**2 * Z1**2 * A)**(1/6)) * (T9)**(5/6)

    E_min = max(0, E0 - 2 * delta_E)
    E_max = E0 + 2 * delta_E

    # Integrate using the sampled S-factor polynomial
    integral, _ = quad(lambda E: np.exp(-dpieta(E)) * Spoly(E, S0, S1, S2, S3) * np.exp(-E/(k_B*T9)), E_min, E_max)
    return prefactor * integral

# Generate samples and calculate rates
S_samples = sample_S_params(100)
rates_samples = []

for S in S_samples:
    rates_i = [Nasv_with_S(T, S) for T in T_range_GK]
    rates_samples.append(rates_i)

rates_samples = np.array(rates_samples)

# Percentile bounds
lower = np.percentile(rates_samples, 16, axis=0)
upper = np.percentile(rates_samples, 84, axis=0)

# Plot with error band
plt.figure(figsize=(15, 7))
plt.plot(T_range_GK, rates, color='darkred', label='Nominal Rate')
plt.fill_between(T_range_GK, lower, upper, color='salmon', alpha=0.4, label='±1σ Range')
plt.xlabel('Temperature (GK)')
plt.ylabel(r'N$_A\langle\sigma v\rangle$ (cm$^3$ mol$^{-1}$ s$^{-1}$)')
plt.title('Reaction Rate vs Temperature for D(p,γ)³He with Uncertainty')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# In[ ]:




