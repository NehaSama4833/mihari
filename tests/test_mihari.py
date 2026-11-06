import numpy as np
from scipy.stats import multivariate_normal
from mihari import policies as P

N, K = 500, 20


def false_rate(pol, k=K, trials=20000, seed=3, p=0.02):
    c1, c2, n = P.simulate(p, p, N, k, trials, np.random.default_rng(seed))
    return pol(c1, c2, n)[0].mean()


def test_fixed_horizon_false_rollback_is_about_alpha():
    assert abs(false_rate(P.fixed_horizon) - 0.05) < 0.012


def test_peeking_inflates_false_rollbacks_well_above_alpha():
    assert false_rate(P.peeking) > 0.15


def test_two_look_gaussian_peeking_matches_exact_bivariate_normal():
    rng = np.random.default_rng(0)
    steps = rng.standard_normal((200000, 2))
    S = steps.cumsum(1) / np.sqrt(np.arange(1, 3))
    sim = (S > P.Z_ONE_SIDED).any(1).mean()
    rho = np.sqrt(0.5)
    exact = 1 - multivariate_normal(mean=[0, 0], cov=[[1, rho], [rho, 1]]).cdf([P.Z_ONE_SIDED] * 2)
    assert abs(sim - exact) < 0.004
    assert 0.08 < exact < 0.09


def test_calibrated_boundaries_hold_alpha_on_fresh_simulation():
    c = P.calibrate(P.constant_boundary, 1.0, 5.0, 0.02, N, K, 30000, seed=1)
    assert abs(false_rate(P.constant_boundary(c), seed=999) - 0.05) < 0.012


def test_msprt_false_rollback_not_above_alpha_plus_noise():
    assert false_rate(P.msprt(0.02)) < 0.05 + 0.01


def test_all_policies_detect_a_large_regression():
    c1, c2, n = P.simulate(0.02, 0.08, N, K, 5000, np.random.default_rng(4))
    for pol in (P.peeking, P.msprt(0.02), P.fixed_horizon, P.sprt_known_baseline(0.02, 0.04)):
        assert pol(c1, c2, n)[0].mean() > 0.95


def test_first_crossing_hand_case():
    hit, first = P.first_crossing(np.array([[False, True, True], [False, False, False]]))
    assert list(hit) == [True, False] and list(first) == [2, 3]
