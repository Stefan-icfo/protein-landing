"""Final numbers and figure data for draft 3: 40 tubes with random non-idealities (batch_random.py),
tension steps 3e-5 per landing, slit window +-0.05 L, 30 landings per tube, every 4th landing a reference."""
import json, numpy as np
from analyze_tension import stats
R = json.load(open('random_rA.json')) + json.load(open('random_rB.json'))
o, e, M = stats(R)
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h
k = int(round(o['mis'] * o['n']))
print(f"mechanical: n={o['n']} mean {o['mean']:+.1f} SD {o['sd']:.1f} Da, misassigned {k}/{o['n']} (95% CI {wilson(k, o['n'])[0]:.3f}-{wilson(k, o['n'])[1]:.3f})")
per_tube = [stats([r])[0]['sd'] for r in R]
print(f"per-tube SD {min(per_tube):.1f}-{max(per_tube):.1f} Da (median {np.median(per_tube):.1f})")
rng = np.random.default_rng(0)
tot = e + rng.normal(0, np.hypot(7, 5), e.size)          # protein-specific bounds, 7 Da and 5 Da
mis = ((M + tot) > 8604.8) != (M > 8604.8); k2 = int(mis.sum())
print(f"with protein-specific terms: SD {tot.std():.1f} Da, misassigned {k2}/{e.size} (95% CI {wilson(k2, e.size)[0]:.3f}-{wilson(k2, e.size)[1]:.3f})")
bins = np.arange(8480, 8740, 10)
x = M + tot
assert ((x >= bins[0]) & (x < bins[-1])).all()
h_ub, _ = np.histogram(x[M < 8600], bins); h_p, _ = np.histogram(x[M > 8600], bins)
json.dump(dict(bins=bins.tolist(), ub=h_ub.tolist(), pub=h_p.tolist(), n_ub=int((M < 8600).sum()), n_pub=int((M > 8600).sum()),
               sd_mech=float(o['sd']), sd_tot=float(tot.std()), mis=float(k2 / e.size), n_tubes=len(R)), open('hist_final.json', 'w'))
print('ub ', h_ub.tolist()); print('pub', h_p.tolist()); print('max bin', int(max(h_ub.max(), h_p.max())))
