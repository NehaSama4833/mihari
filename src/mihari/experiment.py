"""python -m mihari.experiment -> results/report.md"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from scipy.stats import multivariate_normal
from . import policies as P

P0, N_LOOK, K, TRIALS = 0.02, 500, 20, 20000
SCEN = [0.02, 0.03, 0.04, 0.06]


def main():
    out = Path("results"); out.mkdir(exist_ok=True)
    c_pock = P.calibrate(P.constant_boundary, 1.0, 5.0, P0, N_LOOK, K, 40000, seed=101)
    c_obf = P.calibrate(P.obrien_fleming, 1.0, 8.0, P0, N_LOOK, K, 40000, seed=102)
    pols = {"fixed horizon (test once at look 20)": P.fixed_horizon,
            "naive peeking (z>1.645 every look)": P.peeking,
            f"constant boundary z>{c_pock:.2f} (Pocock-style, calibrated)": P.constant_boundary(c_pock),
            f"O'Brien-Fleming-style (c={c_obf:.2f}, calibrated)": P.obrien_fleming(c_obf),
            "mixture SPRT (always-valid, tau=0.02)": P.msprt(0.02),
            "SPRT vs known baseline (H1 = 4%)": P.sprt_known_baseline(P0, 0.04)}
    L = ["# mihari results\n",
         f"Setup: baseline error rate {P0:.0%}, {N_LOOK} requests per arm per look, up to {K} looks ({N_LOOK*K} per arm), "
         f"{TRIALS} simulated rollouts per cell (Monte Carlo SE of a rate near 5% is about 0.2 points). One-sided alpha 0.05. "
         "Constants for the two calibrated rules were fitted on separate simulations under no regression.\n",
         "Cells: probability of rolling back (mean look at which it rolled back).\n",
         "| policy | " + " | ".join(f"canary {p:.0%}" + (" (no regression)" if p == P0 else "") for p in SCEN) + " |",
         "|---|" + "---|" * len(SCEN)]
    rng = np.random.default_rng(7)
    data = {p: P.simulate(P0, p, N_LOOK, K, TRIALS, rng) for p in SCEN}
    for name, pol in pols.items():
        row = []
        for p in SCEN:
            hit, first = pol(*data[p])
            row.append(f"{hit.mean():.1%} ({first[hit].mean():.1f})" if hit.any() else f"{hit.mean():.1%} (-)")
        L.append(f"| {name} | " + " | ".join(row) + " |")
    L.append("")
    L.append("## Peeking inflation vs number of looks (no regression), compared with exact Gaussian theory for two looks\n")
    L.append("| looks | simulated false rollback (binomial) |")
    L.append("|---|---|")
    for k in (1, 2, 5, 10, 20):
        c1, c2, n = P.simulate(P0, P0, N_LOOK, k, 40000, np.random.default_rng(55))
        L.append(f"| {k} | {P.peeking(c1, c2, n)[0].mean():.1%} |")
    rho = np.sqrt(0.5)
    exact2 = 1 - multivariate_normal(mean=[0, 0], cov=[[1, rho], [rho, 1]]).cdf([P.Z_ONE_SIDED, P.Z_ONE_SIDED])
    L.append(f"\nExact two-look value for equally spaced Gaussian looks (corr = sqrt(1/2)): {exact2:.1%}.\n")
    (out / "report.md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
