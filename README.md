# Event-resolved weighing of single landed proteins on a carbon-nanotube resonator

A Monte Carlo of the protocol: proteins land one at a time through a slit onto the centre of a nanotube
resonator at 4 K. Three mode frequencies and a gate sweep are recorded before and after each landing,
and the mass is recovered from the frequency shifts. The question it answers is how well single-molecule
masses can be recovered once realistic frequency noise, nuisance shifts and an imperfectly known tube are
included.

Author: Stefan Forstner (ICFO). Licence: MIT (see `LICENSE`).

## Reproducing the results

Requirements: Python ≥ 3.9, numpy and scipy (`pip install -r requirements.txt`). All scripts run from
this folder and read and write their result files here. The commands below work on Windows, Linux and
macOS (use `python3` instead of `python` where that is the interpreter's name).

```bash
# results table from the stored runs (seconds)
python summarize.py

# main result: 15 Da mechanical, 17 Da with protein-specific terms, 0.8% misassigned; writes hist_final.json
python final_numbers2.py

# re-run the main case (A) from scratch: 40 tubes, about 6-20 s per tube on one core
python batch_random.py rA 600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619
python batch_random.py rB 620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639
python final_numbers2.py
```

The other cases are re-run the same way with the scripts listed under Files. Every run is seeded.
Re-runs reproduce the statistics, but not every event to the last dalton: the position search can
settle on a slightly different solution when numpy, scipy or the LAPACK build differs. With a newer
library stack than the one that produced the stored files, case A re-ran to 15.2 Da mechanical, and 17.1 Da with 0.9% misassigned with the protein-specific terms (stored: 15.0 Da, 17.0 Da, 0.8%; Python 3.13, numpy 2.5, scipy 1.18).

## Model

**Truth (the "real" tube).** A clamped-clamped beam under tension, finite differences with 399 interior
nodes.
- Geometry and material: L = 300 nm, d = 2 nm, graphene areal density, E = 1 TPa, wall 0.34 nm, T0 = 0.5 nN.
- Unloaded modes: 0.84, 2.00 and 3.62 GHz. Tube mass 0.86 MDa.
- One ubiquitin at the centre shifts the modes by −1.19%, −0.002% and −0.98%.
- Each tube departs from the ideal beam in ways the analysis does not know about (`batch_random.py`):
  - a compliant contact at a random end: bending rigidity × U(0.7, 1.0) over 3% of the length;
  - a linear mass and rigidity gradient of U(−1.5%, +1.5%) from end to end;
  - a contamination lump of U(0, 3) kDa at a position drawn from U(0.05, 0.95) L.
- The fixed "mild" tube used in some cases has a contact softened to × 0.8, a 1% gradient and a 2 kDa
  lump at 0.7 L.

**Landings.**
- Ubiquitin (8564.8 Da) or pSer65-ubiquitin (8644.8 Da), 50:50.
- Each protein is modelled as a Gaussian mass footprint, 1.5 nm (1σ).
- Positions are uniform within ±0.05 L of the centre, which corresponds to a 30 nm slit.
- 30 landings per tube, about 30% of the tube's mass.

**Nuisances, every landing.**
- Frequency noise: 2×10⁻⁶ per reading (Allan deviation at 1 s). Each mode is read as a five-point gate
  sweep of ±0.4 V, and the vertex comes from a parabola fit, while the truth curve also has a quartic term.
- Charge offsets: the vertex of the frequency–gate curve moves by N(0, 0.1 V).
- Tension steps: N(0, 3×10⁻⁵) relative change, cumulative.
- Background molecules: one molecule of 18–44 Da lands at a random position in 1% of the steps.

**Analysis.**
- A uniform, symmetric beam. Bending rigidity and tension are fitted to the three unloaded frequencies,
  which are measured with 0.3 ppm noise.
- Loaded frequencies are computed in a 20-mode modal basis.
- Inversion: mass and position from the three relative frequency shifts, with the tension held fixed.
  The position is found by a scan along the whole tube, refined by Gauss–Newton.
- Earlier landings stay in the model at the nominal mass of 8604.8 Da and their estimated positions.

**Per-tube calibration.**
- Every fourth landing (k = 0, 4, …, 28) is a reference of known mass, standing for mass-selected
  CsI cluster ions interleaved with the proteins.
- For each tube, a straight-line fit of the reference errors against landing index removes the tube's
  offset and drift. This offset is −37 Da on average, mostly from the unknown mode shapes.
- Results are evaluated on the remaining 22 landings per tube.

**Protein-specific terms.**
- Added as Gaussian noise of 7 Da (deposited species) and 5 Da (mass scale).
- These are assumed bounds, not simulated physics.

## Results (`python summarize.py`, also in `results_summary.txt`)

"Mech." is the residual standard deviation after the per-tube correction. "+protein" adds the 7 and 5 Da
terms. "Misassigned" uses a threshold at 8604.8 Da.

| Case | Tubes | Molecules | Mech. SD (Da) | Misassigned | +protein SD | Misassigned |
|---|---|---|---|---|---|---|
| **A: random tubes, slit ±0.05 L, tension steps 3×10⁻⁵ (main result)** | 40 | 880 | **15.0** (8–23 per tube) | 0.5% | **17.0** | **0.8%** (95% CI 0.4–1.6%) |
| B: as A, landing window ±0.10 L | 40 | 880 | 31.3 | 9.4% | 32.5 | 11.0% |
| C: "mild" tube, slit, tension steps 3×10⁻⁵ | 40 | 880 | 19.8 | 2.2% | 21.5 | 3.1% |
| D: "mild" tube, slit, no tension steps | 20 | 440 | 18.6 | 1.8% | 20.4 | 1.8% |
| E: "mild" tube, slit, tension steps 1×10⁻⁴ | 20 | 440 | 36.5 | 10.5% | 37.6 | 12.5% |
| F: "mild" tube, tension fitted as a third unknown, no steps | 20 | 440 | 127 | 36% | – | – |
| G: "mild" tube, 13.7 kDa pair (RNase A and +80 Da), 20 landings | 20 | 300 | 32.4 (0.24% of mass) | 10% | – | – |
| H: ideal tube, no noise or nuisances (analysis floor) | 8 | 176 | 8.7 | 0 | – | – |
| I: ideal tube, frequency noise only | 8 | 176 | 9.4 | 0 | – | – |
| J: ideal tube, charge offsets only | 8 | 176 | 8.7 | 0 | – | – |
| K: "mild" tube, no noise or nuisances | 8 | 176 | 16.3 | 0.6% | – | – |

Case L: with landings anywhere along the "mild" tube (80 tubes, `batch.py` without a window), only
10% of the landings are estimated within 0.05 L of the centre. That leaves about 0.35 usable molecules
per tube before accumulated loads spoil the inversion. This is why the design uses a slit.

**Frequency stability** (`noise_sweep.py`, case A conditions, 20 tubes each, seeds 600–619):

| Allan deviation at 1 s | Mech. SD (Da) | +protein SD (Da) | Misassigned (+protein) |
|---|---|---|---|
| 2×10⁻⁶ (case A) | 14.4 | 16.1 | 0.2% |
| 5×10⁻⁶ | 15.8 | 17.3 | 0.9% |
| 1×10⁻⁵ | 17.9 | 19.6 | 1.8% |
| 2×10⁻⁵ | 26.9 | 27.9 | 7.5% |

**Slit width and stability together** (`window_noise.py`, same conditions, mechanical SD):

| Slit | Allan 1×10⁻⁵ | Allan 2×10⁻⁶ | Allan 2×10⁻⁷ |
|---|---|---|---|
| ±0.05 L (30 nm) | 17.9 Da | 14.4 Da | 14.6 Da |
| ±0.025 L (15 nm) | 16.0 Da | 12.2 Da | 12.1 Da |

Ten times better frequency stability than 2×10⁻⁶ does not help: the error is set by the mode shapes, not
the readout. **Fitting precision** (`refine_test.py`, 4 tubes, no noise): a continuous least-squares
refinement of each landing lowers the ideal-tube floor (case H) from about 8 to 4 Da, so part of that
floor is numerical. On the "mild" tube the error stays at about 18 Da either way.

## What this does and does not show

**What the simulation shows:**
- With landings confined to the centre, the protocol recovers single-protein masses to about 15 Da
  (mechanical, 1σ) for ubiquitin on a 0.9 MDa tube. Adding the assumed protein-specific terms gives
  about 17 Da, and fewer than 2% of molecules are misassigned between states 80 Da apart.
- Decomposition (H–K):
  - the analysis itself leaves 4–9 Da;
  - frequency noise adds about 3 Da in quadrature;
  - the gate sweep removes the charge offsets;
  - the unknown mode shapes add the rest.
- The residual grows roughly in proportion to the landed mass (G), so a relative target (≤0.25% of the
  mass) describes it better than a fixed number of daltons.
- Tension must not be a free parameter (F). Tension steps must stay at about 3×10⁻⁵ per landing or
  below; steps of 10⁻⁴ add about 30 Da (E).

**Assumptions to check against real devices:**
- Mode 3 at 3.6 GHz must be readable.
- The Allan deviation is 2×10⁻⁶ at 1 s (up to 1×10⁻⁵ is tolerated), and drift is slow enough for the
  references to follow.
- Charge offsets are about 0.1 V per landing, and tension steps are 3×10⁻⁵ or less.
- The tube's departures from the ideal beam are of the sizes assumed above. Case C, a fixed softer
  contact, gives 20 Da, so the result depends on the tube.

**Not modelled:**
- mode coupling and nonlinearity;
- single-electron softening of the frequency–gate curve (the model curve is smooth);
- diffusion or desorption of landed proteins;
- damping or stiffness changes caused by the protein;
- a non-uniform beam profile across the slit;
- a different landing distribution for the references than for the proteins.

## Files

**Models:**
- `nanotube_balance_mc.py`: truth beam, analysis model, gate sweep, inversions (`invert`,
  `invert_fixed_tension`).
- `calibrated_mc.py`: beam builder with imperfections (`unit_build`, `freqs2`). Its parametric
  calibration with C60 (`run_tube`) was an earlier attempt that did not work and is not used.
- `final_sim.py`: `make_truth('mild' | 'ideal')`. Its `run_tube` (C60 scale calibration with landings
  along the whole tube) is an earlier attempt and is not used.

**Runs:**

| Script | Cases | Usage |
|---|---|---|
| `batch_random.py` | A, B | `python batch_random.py <tag> <seed,seed,...> [sigma_tau] [window]` |
| `batch_tension.py` | C–F | `python batch_tension.py <sigma_tau> fixed\|free <tag> [seeds]` |
| `batch.py` | D (no tension steps), L | `python batch.py mild <seeds> <tag> [window]` |
| `heavy.py` | G | `python heavy.py <seeds>` |
| `decompose.py` | H–K | `python decompose.py ideal_clean\|noise_only\|charge_only\|tube_only` |
| `noise_sweep.py` | stability sweep | `python noise_sweep.py <allan_dev> <tag> <seeds>` |
| `window_noise.py` | slit × stability | `python window_noise.py <allan_dev> <window> <tag>` |
| `refine_test.py` | fitting precision | `python refine_test.py ideal mild` |

**Analysis:**
- `analyze_tension.py`: per-tube reference correction and statistics, used by all summaries.
- `summarize.py`: prints the results table.
- `final_numbers2.py`: main-result numbers and histogram data (`hist_final.json`).
- `final_numbers.py`: earlier analysis of case D (`hist.json`).

**Stored results:**
- Case A: `random_rA.json`, `random_rB.json`. Case B: `random_w10A.json`, `random_w10B.json`.
- Cases C–F: `tension_*.json`. D and L: `batch_win05.json`, `batch_win10.json`, `batch_mildA.json`,
  `batch_mildB.json`. G: `batch_heavy.json`.
- Cases H–K: `decomp_*.json`.
- Stability sweep: `noise_*.json`. Slit × stability: `wn_*.json`.
