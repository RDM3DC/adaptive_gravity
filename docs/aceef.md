# ACEEF cosmology extension

This module introduces optional Adaptive Conductance Extended Energy Framework (ACEEF) pieces into the Adaptive Gravity codebase.

## Equations

Galaxy runs keep the original adaptive conductance field:

\[\dot G_c(r,t)=\alpha\,\rho(r)^2-\mu\,G_c(r,t)+D\,\partial_{rr}G_c\]

Cosmology adds a slow hyperspatial channel and a holographic vacuum term:

\[\dot G_h = \beta_h\,(\rho_\gamma H)\,e^{-a/a_c}-\gamma_h G_h\]
\[\dot\rho_\Lambda = \kappa\,\frac{d}{dt}\Big(\frac{a^2}{L_p^2}\rho_{\rm holo}\Big)-2\kappa\,\frac{L_p^2}{a^3}\,\dot a\]

The background expansion uses the combined field \(G_{\rm eff}=G_c+\lambda_h G_h\):

\[H^2=\frac{8\pi}{3}G_{\rm eff}(\rho_m+\rho_r)+\frac{\Lambda}{3}\,,\qquad \Lambda=8\pi\rho_\Lambda.\]

## Switches

- `lambda_h` – weight of the hyperspatial field. Defaults to 0 to preserve galaxy fits.
- `beta_h`, `gamma_h`, `a_c` – govern the evolution of `G_h`.
- `kappa` – holographic coupling; 0 turns the feature off.

## Baseline run

A minimal parameter file lives at `cosmo/params/baseline.json`. It turns off all new couplings and reproduces the prior constant‑G behaviour.

Run it with:

```bash
python -m adaptive_gravity.cosmology cosmo/params/baseline.json runs/cosmo_baseline/
```

This writes an `evolution.csv` table to the output directory.

## Tested vs speculative

- `G_h` and `\rho_\Lambda(t)` are implemented and can be enabled via the parameters above.
- Entanglement amplifiers and chaotic terms from ACEEF are **not** yet coded; the structure allows adding them later.

The baseline command above reproduces the nominal evolution and serves as a starting point for further experimentation.
