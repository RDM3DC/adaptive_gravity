"""Minimal cosmology module with optional ACEEF terms."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import numpy as np


def _compute_H(a: float, Gc: float, Gh: float, rhoL: float, pars: Dict[str, float]) -> float:
    """Return the Hubble rate given state variables and parameters."""
    H2 = (8 * np.pi / 3.0) * (Gc + pars["lambda_h"] * Gh) * (
        pars["Om"] / a**3 + pars["Or"] / a**4
    )
    H2 += (8 * np.pi * rhoL) / 3.0
    return float(np.sqrt(max(H2, 0.0)))


def rhs(t: float, y: np.ndarray, pars: Dict[str, float]) -> np.ndarray:
    """Time derivatives for [a, Gc, Gh, rhoL]."""
    a, Gc, Gh, rhoL = y
    H = _compute_H(a, Gc, Gh, rhoL, pars)
    rho_m = pars["Om"] / a**3
    rho_g = pars["Or"] / a**4
    dGc = pars["alpha_c"] * rho_m**2 - pars["mu_c"] * Gc
    dGh = pars["beta_h"] * (rho_g * H) * np.exp(-a / pars["a_c"]) - pars["gamma_h"] * Gh
    dRhoL = (
        pars["kappa"] * ((2 * a) * pars["rho_holo"] / pars["Lp2"] * H)
        - 2 * pars["kappa"] * pars["Lp2"] * pars["rho_holo"] * H / a**3
    )
    return np.array([a * H, dGc, dGh, dRhoL])


def run(params_path: str, out_dir: str = "runs/cosmo") -> str:
    """Integrate the cosmological system and write a CSV of the evolution."""
    cfg = json.loads(Path(params_path).read_text())
    pars = cfg
    init = cfg.get("init", {"a": 1.0, "Gc": 1.0, "Gh": 0.0, "rhoL": 0.0})
    tcfg = cfg.get("t", {"t_start": 0.0, "t_end": 1.0, "n_steps": 100})

    y = np.array([init["a"], init["Gc"], init["Gh"], init["rhoL"]], dtype=float)
    t = float(tcfg.get("t_start", 0.0))
    dt = (tcfg["t_end"] - tcfg["t_start"]) / tcfg["n_steps"]

    data = []
    for _ in range(tcfg["n_steps"]):
        H = _compute_H(y[0], y[1], y[2], y[3], pars)
        data.append([t, y[0], y[1], y[2], y[3], H])
        k1 = rhs(t, y, pars)
        k2 = rhs(t + 0.5 * dt, y + 0.5 * dt * k1, pars)
        k3 = rhs(t + 0.5 * dt, y + 0.5 * dt * k2, pars)
        k4 = rhs(t + dt, y + dt * k3, pars)
        y += (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
    H = _compute_H(y[0], y[1], y[2], y[3], pars)
    data.append([t, y[0], y[1], y[2], y[3], H])

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    import pandas as pd

    df = pd.DataFrame(data, columns=["t", "a", "Gc", "Gh", "rhoL", "H"])
    csv_path = out / "evolution.csv"
    df.to_csv(csv_path, index=False)
    return str(csv_path)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        raise SystemExit("Usage: python -m adaptive_gravity.cosmology <params.json> [out_dir]")
    params = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "runs/cosmo"
    print(run(params, out_dir))
