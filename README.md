# Adaptive Gravity (ARP-Driven) – Galaxy Rotation Simulation

This project implements an **adaptive gravity model** inspired by the Adaptive Resistance Principle (ARP), applied to galactic rotation curves.

## 📐 Model Summary

- Gravity evolves dynamically as:
  \[ \dot{G}(r,t) = \alpha \rho(r)^2 - \mu G(r,t) \]
- Enforces a pivot at radius \( R_0 = 8 \text{ kpc} \) where all rotation curves cross (set by equilibrium \( G = G_N \)).
- Memory decay \( \mu = 1/T_0 \), where \( T_0 = 2\pi R_0 / V_0 \).
- Uses realistic Milky Way-style mass model:
  - **Hernquist bulge**
  - **Exponential disk** (with vertical scale height for 3D density)

## 📊 Outputs

- `bulge_disk_adaptive_rotation.csv` — Full output table
- `rotation_curve_bulge_disk.png` — Rotation curve at \( t = 2T_0 \)
- `acceleration_profile_adaptive.png` — Acceleration curve with MOND \( a_0 \) overlay

## 🧠 How to Use

1. Run `adaptive_gravity_sim.py` to reproduce all outputs.
2. Edit `μ`, `D`, or mass model parameters to explore different behaviors.
3. Compare with observed galaxies (e.g., SPARC) for fitting tests.

## 📎 Files

- `adaptive_gravity_sim.py` — Python simulation code
- `rotation_curve_bulge_disk.png` — Output figure (velocity)
- `acceleration_profile_adaptive.png` — Output figure (acceleration)
- `bulge_disk_adaptive_rotation.csv` — Results table
- `README.md` — This file

---

Built by @RDM3DC • Aug 2025
