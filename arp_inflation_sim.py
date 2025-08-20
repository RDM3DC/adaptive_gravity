# arp_inflation_sim.py
import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11         # Gravitational constant
rho_vac_0 = 3.09e103    # Unmodified quantum vacuum energy density (J/m^3)
c = 299792458           # Speed of light (m/s)

# ARP Parameters
mu = 1e28
alpha = 1e-38           # Controls decay of vacuum energy

def rho_vac(t):
    """Exponential decay of vacuum energy over time"""
    return rho_vac_0 * np.exp(-alpha * t)

def G_eff(t):
    """Adaptive gravitational constant based on energy density"""
    rho = rho_vac(t)
    return G * np.exp(-mu / (rho + 1e-30))

def friedmann_rhs(t, a):
    """Right-hand side of Friedmann equation for a(t)"""
    return np.sqrt((8 * np.pi * G_eff(t) * rho_vac(t)) / 3) * a

# Time and integration setup
t_max = 1e-35          # End time in seconds
dt = 1e-39             # Time step
t_vals = [0.0]
a_vals = [1e-30]       # Start with a tiny universe

# Runge-Kutta 4th order integration
while t_vals[-1] < t_max:
    t = t_vals[-1]
    a = a_vals[-1]

    k1 = dt * friedmann_rhs(t, a)
    k2 = dt * friedmann_rhs(t + dt/2, a + k1/2)
    k3 = dt * friedmann_rhs(t + dt/2, a + k2/2)
    k4 = dt * friedmann_rhs(t + dt, a + k3)
    
    a_next = a + (k1 + 2*k2 + 2*k3 + k4) / 6
    t_next = t + dt

    t_vals.append(t_next)
    a_vals.append(a_next)

# Convert to numpy arrays for plotting
t_vals = np.array(t_vals)
a_vals = np.array(a_vals)
rho_vals = rho_vac(t_vals)
G_vals = G_eff(t_vals)

# Plot results
plt.figure(figsize=(12,6))
plt.subplot(1,3,1)
plt.plot(t_vals, a_vals)
plt.title("Scale Factor a(t)")
plt.xlabel("Time (s)")
plt.ylabel("a(t)")

plt.subplot(1,3,2)
plt.plot(t_vals, rho_vals)
plt.title("Vacuum Energy Density")
plt.xlabel("Time (s)")
plt.ylabel(r"$\rho_{vac}(t)$ (J/m$^3$)")

plt.subplot(1,3,3)
plt.plot(t_vals, G_vals)
plt.title("Effective Gravity $G_{eff}(t)$")
plt.xlabel("Time (s)")
plt.ylabel(r"$G_{eff}(t)$")

plt.tight_layout()
plt.savefig("arp_inflation_plots.png")
plt.show()
