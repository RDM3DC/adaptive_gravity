# quantum_foam_arp.py
import numpy as np
import matplotlib.pyplot as plt

# Constants
hbar = 1.0545718e-34  # Js
c = 299792458         # m/s
G = 6.67430e-11        # m^3 kg^-1 s^-2
l_p = np.sqrt(hbar * G / c**3)
t_p = np.sqrt(hbar * G / c**5)

# ARP Parameters
mu = 1e20
alpha = 1e-40

# Simulation domain
N = 10000  # Number of Planck cells
L = N * l_p  # Total simulated length (1D foam)
x = np.linspace(0, L, N)

# Random curvature fluctuation field (1/r^2-like bumps with noise)
np.random.seed(42)
curvature = np.abs(np.random.normal(loc=1/l_p**2, scale=1e60, size=N))

# Adaptive G field
def G_eff(curv):
    return G * np.exp(-mu * curv)

G_field = G_eff(curvature)

delta_E = hbar / t_p  # Energy fluctuation per cell (uncertainty limit)

# Effective energy density from vacuum
energy_density = (G - G_field) * delta_E / (l_p**3)

# Average vacuum energy density
rho_vac = np.mean(energy_density)

print(f"Estimated vacuum energy density from ARP foam: {rho_vac:.3e} J/m^3")

# Compare with observed dark energy density (Lambda)
rho_lambda = 6.91e-10  # J/m^3
print(f"Observed dark energy density: {rho_lambda:.3e} J/m^3")
print(f"Ratio (ARP / Lambda): {rho_vac / rho_lambda:.3e}")

# Plot distribution
plt.figure(figsize=(8,5))
plt.hist(energy_density, bins=100, color='purple', alpha=0.7, log=True)
plt.axvline(rho_lambda, color='red', linestyle='--', label='Observed Dark Energy')
plt.axvline(rho_vac, color='green', linestyle='--', label='ARP Vacuum Mean')
plt.xlabel("Energy Density (J/m^3)")
plt.ylabel("Frequency (log scale)")
plt.title("ARP Quantum Foam Energy Density Fluctuations")
plt.legend()
plt.tight_layout()
plt.savefig("arp_quantum_foam_density.png")
plt.show()
