import numpy as np
import matplotlib.pyplot as plt

# Constants and settings
G = 4.30091e-6  # kpc·(km/s)^2 / Msun
R0_kpc = 8.0
V0_kms = 220.0
T0 = 2 * np.pi * R0_kpc / V0_kms
mu = 1 / T0

# Bulge (Hernquist)
Mb = 0.9e10  # Msun
a_b = 0.7  # kpc

def M_bulge(r): return Mb * r**2 / (r + a_b)**2
def rho_bulge(r): return (Mb / (2 * np.pi)) * a_b / (r * (r + a_b)**3 + 1e-8)

# Disk (Exponential)
Md = 5e10  # Msun
Rd = 2.6  # kpc
z0 = 0.3  # scale height

def Sigma_d(R): return (Md / (2 * np.pi * Rd**2)) * np.exp(-R / Rd)
def M_disk(R): return Md * (1 - (1 + R / Rd) * np.exp(-R / Rd))
def rho_disk(R): return Sigma_d(R) / (2 * z0)

# Grid
r = np.linspace(0.01, 40, 1000)
rho = rho_bulge(r) + rho_disk(r)
M_enc = M_bulge(r) + M_disk(r)
rho_R0 = np.interp(R0_kpc, r, rho)
alpha = mu * G / (rho_R0**2)

# Evolve adaptive G
Geff = np.full_like(r, G)
dt = T0 / 500
for _ in range(int(2 * T0 / dt)):
    dG = alpha * rho**2 - mu * Geff
    Geff += dt * dG

v = np.sqrt(np.clip(Geff * M_enc / r, 0, None))
a_kms2_kpc = v**2 / r
a_m_s2 = a_kms2_kpc * 3.24078e-14

# Plot results
plt.plot(r, v)
plt.xlabel("r [kpc]")
plt.ylabel("v [km/s]")
plt.title("Adaptive Gravity Rotation Curve")
plt.savefig("rotation_curve_bulge_disk.png")
plt.close()

plt.plot(r, a_m_s2)
plt.axhline(1.2e-10, color='red', linestyle='--', label='a₀')
plt.xlabel("r [kpc]")
plt.ylabel("a(r) [m/s²]")
plt.title("Acceleration Profile")
plt.legend()
plt.savefig("acceleration_profile_adaptive.png")
plt.close()
