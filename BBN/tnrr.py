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

# # integrand  
#   integrand <- function(E, T9){
#                exp(-dpieta(E)) * Spoly(E, alpha, beta,delta,gamma) * 
#                			   exp(-E/(0.086173324*T9))
#                              }
def integrand(E, T9):
    """
    integrand for the integral
    """
    return np.exp(-dpieta(E)) * Spoly(E, S0, S1, S2, S3) * np.exp(-E/(k_B*T9))

  # reaction rate
#   Nasv <- function(Temp){
#          sqrt(8/ (pi * mue * muc2)) * clight * avog * 1e-24 *
#          ((0.086173324 * Temp)^{-3/2}) * 
#          integrate(integrand, lower = 1e-6,
# 						upper = 1.9,abs.tol = 0L,
# #						upper = 2,abs.tol = 0L,
# 						T9 = Temp)$value
def Nasv(T):
    """
    Reaction rate
    """
    sqrt_term = np.sqrt(8 / (np.pi * mue * muc2))
    avog = N_A
    prefactor = sqrt_term * clight * avog * 1e-24 *  ((k_B * T)**(-3/2))
    integral, _ = quad(integrand, 1e-6, 2, args=(T,))
    return prefactor * integral

print("THermp",Nasv(.1))

def e0(Z1,Z2,A,T9):
    return 0.1220*((Z1**2)*(Z2**2)*A)**(1/3) * T9**(2/3) 

def dele0(Z1,Z2,A,T9):
    return 0.2368*((Z1**2)*(Z2**2)*A)**(1/6) * T9**(2/6) 

print(e0(Z0,Z1,mue,10))

print(dele0(Z0,Z1,mue,10))