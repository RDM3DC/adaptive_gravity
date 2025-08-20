# arp_reheating_sim.py
import numpy as np
import matplotlib.pyplot as plt

# Constants
c = 299792458          # m/s
hbar = 1.0545718e-34   # Js
G = 6.67430e-11        # m^3 kg^-1 s^-2
kB = 1.380649e-23      # J/K

# ARP parameters
mu = 1e28      # decay suppression
alpha = 1e-38  # vacuum energy inertia

# Time range: from end of inflation to early Big Bang (~10^-32 to 10^-10 s)
t = np.logspace(-32, -10, 500)

# Initial conditions
rho_vac0 = 3.09e103   # initial vacuum energy (J/m³)
T_rad0 = 1e9          # starting radiation temp (K)

# Decay function for vacuum energy -> radiation
def rho_vac(t):
    decay = np.exp(-mu * t) * np.cos(1e32 * t)**2  # oscillating decay
    return rho_vac0 * decay

def rho_rad(t):
    return rho_vac0 - rho_vac(t)  # energy transferred from vacuum to radiation

def temperature_rad(rho_r):
    """Estimate radiation temperature from energy density"""
    a_rad = (np.pi**2 / 15) * (kB**4) / (hbar**3 * c**5)
    return (rho_r / a_rad)**0.25

# Compute values
rho_v = rho_vac(t)
rho_r = rho_rad(t)
T_r = temperature_rad(rho_r)

# Plot energy densities
plt.figure(figsize=(10,5))
plt.loglog(t, rho_v, label="Vacuum Energy (ρ_vac)")
plt.loglog(t, rho_r, label="Radiation Energy (ρ_rad)")
plt.xlabel("Time (s)")
plt.ylabel("Energy Density (J/m³)")
plt.title("ARP Reheating: Energy Transfer After Inflation")
plt.legend()
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.savefig("arp_reheating_energy.png")
plt.show()

# Plot radiation temperature
plt.figure(figsize=(8,5))
plt.loglog(t, T_r, color="crimson", label="Radiation Temperature")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("ARP Reheating: Radiation Temperature Rise")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.savefig("arp_reheating_temp.png")
plt.show()
