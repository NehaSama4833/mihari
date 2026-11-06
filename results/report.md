# mihari results

Setup: baseline error rate 2%, 500 requests per arm per look, up to 20 looks (10000 per arm), 20000 simulated rollouts per cell (Monte Carlo SE of a rate near 5% is about 0.2 points). One-sided alpha 0.05. Constants for the two calibrated rules were fitted on separate simulations under no regression.

Cells: probability of rolling back (mean look at which it rolled back).

| policy | canary 2% (no regression) | canary 3% | canary 4% | canary 6% |
|---|---|---|---|---|
| fixed horizon (test once at look 20) | 5.1% (20.0) | 99.8% (20.0) | 100.0% (20.0) | 100.0% (20.0) |
| naive peeking (z>1.645 every look) | 21.5% (6.0) | 99.9% (3.5) | 100.0% (1.6) | 100.0% (1.0) |
| constant boundary z>2.39 (Pocock-style, calibrated) | 5.1% (6.9) | 99.2% (6.1) | 100.0% (2.4) | 100.0% (1.2) |
| O'Brien-Fleming-style (c=1.84, calibrated) | 5.2% (15.3) | 99.8% (8.9) | 100.0% (5.1) | 100.0% (3.1) |
| mixture SPRT (always-valid, tau=0.02) | 0.6% (6.5) | 93.1% (9.1) | 100.0% (3.4) | 100.0% (1.5) |
| SPRT vs known baseline (H1 = 4%) | 0.9% (1.5) | 85.3% (5.2) | 100.0% (1.6) | 100.0% (1.0) |

## Peeking inflation vs number of looks (no regression), compared with exact Gaussian theory for two looks

| looks | simulated false rollback (binomial) |
|---|---|
| 1 | 4.9% |
| 2 | 7.9% |
| 5 | 12.9% |
| 10 | 17.5% |
| 20 | 21.7% |

Exact two-look value for equally spaced Gaussian looks (corr = sqrt(1/2)): 8.0%.
