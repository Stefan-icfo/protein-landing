"""Prints the table in README.md from the stored result files."""
import json, numpy as np
from analyze_tension import stats
def load(*fs):
    R = []
    for f in fs: R += json.load(open(f))
    return R
def row(name, R, nominal=8604.8, protein_terms=True):
    o, e, M = stats(R, nominal)
    s = f"{name:<58s} {len(R):>3d} {o['n']:>4d} {o['sd']:>6.1f} {100*o['mis']:>6.1f}"
    if protein_terms:
        rng = np.random.default_rng(0); tot = e + rng.normal(0, np.hypot(7, 5), e.size)
        mis = np.mean(((M + tot) > nominal) != (M > nominal))
        s += f" {tot.std():>7.1f} {100*mis:>6.1f}"
    print(s)
print(f"{'case':<58s} {'tub':>3s} {'mol':>4s} {'SD/Da':>6s} {'mis%':>6s} {'+prot':>7s} {'mis%':>6s}")
row('A  random tubes, slit +-0.05 L, tension steps 3e-5 [MAIN]', load('random_rA.json', 'random_rB.json'))
row('B  random tubes, window +-0.10 L, tension steps 3e-5', load('random_w10A.json', 'random_w10B.json'))
row('C  "mild" tube, slit, tension steps 3e-5', load('tension_t03fixed.json', 'tension_t03fixedB.json', 'tension_t03fixedC.json'))
row('D  "mild" tube, slit, no tension steps', load('batch_win05.json'))
row('E  "mild" tube, slit, tension steps 1e-4', load('tension_t1fixed.json'))
row('F  "mild" tube, slit, tension fitted as 3rd unknown, no steps', load('tension_t0free.json'))
row('G  "mild" tube, slit, 13.7 kDa pair (RNase A, +80 Da)', load('batch_heavy.json'), nominal=13722.0, protein_terms=False)
row('H  ideal tube, no noise, no nuisances (analysis floor)', load('decomp_ideal_clean.json'), protein_terms=False)
row('I  ideal tube, frequency noise only', load('decomp_noise_only.json'), protein_terms=False)
row('J  ideal tube, charge offsets only', load('decomp_charge_only.json'), protein_terms=False)
row('K  "mild" tube, no noise, no nuisances', load('decomp_tube_only.json'), protein_terms=False)
R2 = load('batch_mildA.json', 'batch_mildB.json')
ev = [e for r in R2 for e in r['events']]
xe = np.array([e['xe'] for e in ev]); kk = np.array([e['k'] for e in ev])
print(f"L  'mild' tube, landings along the whole tube ({len(R2)} tubes): fraction estimated within 0.05 L of centre "
      f"{np.mean(np.abs(xe-0.5)<=0.05):.2f}; usable molecules per tube (central, first 5 landings) {np.sum((np.abs(xe-0.5)<=0.05)&(kk<5))/len(R2):.2f}")
per = [stats([r])[0]['sd'] for r in load('random_rA.json', 'random_rB.json')]
print(f"Case A per-tube SD: {min(per):.1f}-{max(per):.1f} Da, median {np.median(per):.1f} Da")
