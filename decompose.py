"""Error decomposition: frequency noise alone (ideal tube, no other nuisance) and tube non-ideality alone (no noise)."""
import json, numpy as np, sys
import nanotube_balance_mc as s, calibrated_mc as c, final_sim as fsim
from analyze_tension import stats
Da, L = s.Da, s.L
UB, PUB, NOM = 8564.8, 8644.8, 8604.8
def tube(seed, kind, sigma_y, sigma_V, p_bg, n_land=30, window=0.05):
    rng = np.random.default_rng(seed)
    Bt = fsim.make_truth(kind)
    f0 = c.freqs2(Bt, c.EI0, c.T0, Bt['m'])
    A = s.Analysis(f0 * (1 + (0.3e-6 if sigma_y > 0 else 0.0) * rng.standard_normal(3)))
    tab = s.sens_table(A)
    loads = []; Mb = A.I.copy(); V0 = 0.0
    fb, Vop = s.sweep(f0, V0, 0.0, sigma_y, rng)
    out = []
    for k in range(n_land):
        M = UB if rng.random() < 0.5 else PUB
        x = L * (0.5 + rng.uniform(-window, window)); loads.append((M * Da, x))
        V0 += sigma_V * rng.standard_normal()
        if rng.random() < p_bg: loads.append((rng.uniform(18, 44) * Da, rng.uniform(0, L)))
        m = Bt['m'].copy()
        for dm, xx in loads: m += dm * s.footprint(Bt['xn'], xx)
        fa, Vop = s.sweep(c.freqs2(Bt, c.EI0, c.T0, m), V0, Vop, sigma_y, rng)
        r = np.log(fa / fb); fb = fa
        me, xe, te, chi = s.invert_fixed_tension(A, Mb, 0.0, r, tab)
        Mb = Mb + A.Mload(NOM * Da, xe * L)
        out.append(dict(k=k, M=M, est=float(me), xt=float(x / L), xe=float(xe), chi=float(chi)))
    return dict(seed=seed, events=out)
if __name__ == '__main__':
    case = sys.argv[1]
    cfg = {'ideal_clean': ('ideal', 1e-9, 0.0, 0.0), 'noise_only': ('ideal', 2e-6, 0.0, 0.0), 'charge_only': ('ideal', 1e-9, 0.1, 0.0), 'tube_only': ('mild', 1e-9, 0.0, 0.0)}[case]
    R = [tube(sd, *cfg) for sd in range(400, 408)]
    json.dump(R, open(f'decomp_{case}.json', 'w'))
    o, e, M = stats(R)
    print(case, {k: round(v, 2) for k, v in o.items()})
