"""Canary decision rules. Control and canary each get n_per_look requests per look; the rule decides at
each look whether to roll back (canary error rate worse) or keep going. Vectorised over trials."""
from __future__ import annotations
import numpy as np
from scipy.stats import norm

ALPHA = 0.05
Z_ONE_SIDED = norm.ppf(1 - ALPHA)


def simulate(p_ctrl, p_can, n_per_look, K, trials, rng):
    c_ctrl = rng.binomial(n_per_look, p_ctrl, size=(trials, K)).cumsum(1)
    c_can = rng.binomial(n_per_look, p_can, size=(trials, K)).cumsum(1)
    n = n_per_look * np.arange(1, K + 1)
    return c_ctrl, c_can, n


def zstat(c_ctrl, c_can, n):
    pc, pk = c_ctrl / n, c_can / n
    pbar = (c_ctrl + c_can) / (2 * n)
    se = np.sqrt(np.maximum(2 * pbar * (1 - pbar) / n, 1e-12))
    return (pk - pc) / se


def first_crossing(crossed):
    """crossed: bool (trials, K). Returns (rolled_back bool, look index 1..K or K if never)."""
    any_ = crossed.any(1)
    first = np.where(any_, crossed.argmax(1) + 1, crossed.shape[1])
    return any_, first


def fixed_horizon(c1, c2, n):
    z = zstat(c1, c2, n)
    cross = np.zeros_like(z, dtype=bool)
    cross[:, -1] = z[:, -1] > Z_ONE_SIDED
    return first_crossing(cross)


def peeking(c1, c2, n):
    return first_crossing(zstat(c1, c2, n) > Z_ONE_SIDED)


def constant_boundary(c):
    return lambda c1, c2, n: first_crossing(zstat(c1, c2, n) > c)


def obrien_fleming(c):
    def f(c1, c2, n):
        K = len(n)
        b = c * np.sqrt(K / np.arange(1, K + 1))
        return first_crossing(zstat(c1, c2, n) > b)
    return f


def msprt(tau=0.02, alpha=ALPHA):
    """Always-valid mixture SPRT for a difference of means with a N(0, tau^2) mixing prior
    (Johari et al. style): roll back when the likelihood ratio >= 1/alpha and canary is worse."""
    def f(c1, c2, n):
        pc, pk = c1 / n, c2 / n
        theta = pk - pc
        pbar = (c1 + c2) / (2 * n)
        V = np.maximum(2 * pbar * (1 - pbar) / n, 1e-12)
        lam = np.sqrt(V / (V + tau ** 2)) * np.exp(tau ** 2 * theta ** 2 / (2 * V * (V + tau ** 2)))
        return first_crossing((lam >= 1 / alpha) & (theta > 0))
    return f


def sprt_known_baseline(p0, p1, alpha=ALPHA):
    """Wald SPRT on the canary alone against a KNOWN baseline rate p0 (H1: rate p1). Needs a trusted p0."""
    a = np.log(1 / alpha)
    def f(c1, c2, n):
        llr = c2 * np.log(p1 / p0) + (n - c2) * np.log((1 - p1) / (1 - p0))
        return first_crossing(llr >= a)
    return f


def calibrate(make_policy, lo, hi, p, n_per_look, K, trials, seed, target=ALPHA):
    """Bisect the policy constant so false rollback under H0 equals target (on its own seed)."""
    rng = np.random.default_rng(seed)
    c1, c2, n = simulate(p, p, n_per_look, K, trials, rng)
    for _ in range(30):
        mid = (lo + hi) / 2
        fa = make_policy(mid)(c1, c2, n)[0].mean()
        lo, hi = (mid, hi) if fa > target else (lo, mid)
    return (lo + hi) / 2
