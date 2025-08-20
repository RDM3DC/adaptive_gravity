from __future__ import annotations
import json, math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

G = 4.30091e-6  # kpc (km/s)^2 / Msun

@dataclass
class MassParams:
    Mb: float = 0.0
    a_bulge: float = 1.0
    Md: float = 0.0
    Rd: float = 2.5
    z0: float = 0.3
    Mg: float = 0.0
    Rg: float = 6.0
    gas_bump_r0: Optional[float] = None
    gas_bump_sigma: float = 1.0
    gas_bump_amp: float = 0.0

@dataclass
class ARPParams:
    alpha: float
    mu: float
    D: float = 0.0
    R0: float = 8.0
    t_units_T0: float = 2.0

def _hernquist_M(r, Mb, a):
    return Mb * (r**2) / (r + a)**2

def _expdisk_Sigma(R, Md, Rd):
    return (Md / (2 * math.pi * Rd**2)) * np.exp(-R / Rd)

def _expdisk_M(R, Md, Rd):
    return Md * (1 - (1 + R / Rd) * np.exp(-R / Rd))

def _gas_density(R, Mg, Rg, bump_r0, bump_sigma, bump_amp, z0):
    Sigma = (Mg / (2 * math.pi * Rg**2)) * np.exp(-R / Rg)
    if bump_r0 is not None and bump_amp != 0.0:
        Sigma = Sigma * (1.0 + bump_amp * np.exp(-0.5 * ((R - bump_r0) / bump_sigma) ** 2))
    rho = Sigma / (2 * z0)
    return rho

def _laplacian_neumann(Garr, dr):
    L = np.zeros_like(Garr)
    L[1:-1] = (Garr[2:] - 2*Garr[1:-1] + Garr[:-2]) / (dr**2)
    L[0] = (Garr[1] - Garr[0]) / (dr**2) * 2.0
    L[-1] = (Garr[-2] - Garr[-1]) / (dr**2) * 2.0
    return L

def run(params_path: str, out_dir: str = "runs/out"):
    cfg = json.loads(Path(params_path).read_text())
    m = MassParams(**cfg["mass"])
    a = ARPParams(**cfg["arp"])
    grid = cfg.get("grid", {"r_min": 0.01, "r_max": 40.0, "nr": 1200})
    anchor = cfg.get("anchor", {"R0": a.R0, "V0": 220.0})

    r = np.linspace(grid["r_min"], grid["r_max"], grid["nr"])
    dr = r[1] - r[0]

    rho_b = (m.Mb / (2 * math.pi)) * m.a_bulge / (r * (r + m.a_bulge) ** 3 + 1e-12)
    rho_d = _expdisk_Sigma(r, m.Md, m.Rd) / (2 * m.z0)
    rho_g = _gas_density(r, m.Mg, m.Rg, m.gas_bump_r0, m.gas_bump_sigma, m.gas_bump_amp, m.z0)
    rho = rho_b + rho_d + rho_g

    M_enc = _hernquist_M(r, m.Mb, m.a_bulge) + _expdisk_M(r, m.Md, m.Rd)
    M_enc += m.Mg * (1 - (1 + r / m.Rg) * np.exp(-r / m.Rg))

    T0 = 2 * math.pi * anchor["R0"] / anchor["V0"]
    rho_R0 = np.interp(a.R0, r, rho)
    alpha = a.mu * G / (rho_R0 ** 2)

    Geff = np.full_like(r, G)
    dt = T0 / 800.0
    steps = int(a.t_units_T0 * T0 / dt)

    def rhs(Garr):
        term = alpha * (rho ** 2) - a.mu * Garr
        if a.D != 0.0:
            term += a.D * _laplacian_neumann(Garr, dr)
        return term

    for _ in range(steps):
        k1 = rhs(Geff)
        k2 = rhs(Geff + 0.5 * dt * k1)
        k3 = rhs(Geff + 0.5 * dt * k2)
        k4 = rhs(Geff + dt * k3)
        Geff += (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

    v = np.sqrt(np.clip(Geff * M_enc / r, 0.0, None))
    a_ms2 = (v**2 / r) * 3.24078e-14

    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    df = pd.DataFrame({
        "r_kpc": r,
        "rho_bulge": rho_b,
        "rho_disk": rho_d,
        "rho_gas": rho_g,
        "rho_total": rho,
        "M_enc_Msun": M_enc,
        "G_eff": Geff,
        "v_kms": v,
        "a_m_s2": a_ms2
    })
    df.to_csv(out / "profile.csv", index=False)

    import matplotlib.pyplot as plt
    plt.figure(figsize=(7,5))
    plt.plot(r, v, label="ARP")
    plt.axvline(a.R0, ls="--", c="gray", lw=1, label=f"pivot R0={a.R0} kpc")
    plt.xlabel("r [kpc]"); plt.ylabel("v [km/s]")
    plt.title("Adaptive Gravity Rotation Curve")
    plt.legend(); plt.tight_layout()
    plt.savefig(out / "rotation_curve.png", dpi=140); plt.close()

    return str(out / "profile.csv"), str(out / "rotation_curve.png")
