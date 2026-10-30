"""Like batch_tension.py (fixed-tension inversion, tension steps sigma_tau per landing), but every tube
has its own random non-ideality, unknown to the analysis:
  compliant contact at a random end: bending rigidity x U(0.7, 1.0) over 3 % of the length
  linear mass/rigidity gradient U(-1.5 %, +1.5 %) end to end
  contamination lump U(0, 3) kDa at a random position in U(0.05, 0.95) L"""
import numpy as np, json, sys, time
import nanotube_balance_mc as s, calibrated_mc as c

Da, L = s.Da, s.L
UB, PUB, NOM = 8564.8, 8644.8, 8604.8

def random_truth(rng):
    k = rng.uniform(0.7, 1.0); g = rng.uniform(-0.015, 0.015)
    if rng.random() < 0.5:
        B = c.unit_build(1.0, k, 1.0, g, softL=0.03, softR=0.0)
    else:
        B = c.unit_build(1.0, 1.0, k, g, softL=0.0, softR=0.03)
    lump, xl = rng.uniform(0, 3000), rng.uniform(0.05, 0.95)
    w = np.exp(-0.5 * ((B['xn'] - xl * L) / 3e-9)**2)
    B['m'] = B['m'] + lump * Da * w / w.sum()
    return B, dict(k=float(k), g=float(g), lump=float(lump), xl=float(xl))

def tube(seed, sigma_tau=3e-5, n_land=30, sigma_y=2e-6, sigma_V=0.1, p_bg=0.01, window=0.05):
    rng = np.random.default_rng(seed)
    Bt, par = random_truth(rng)
    f0 = c.freqs2(Bt, c.EI0, c.T0, Bt['m'])
    A = s.Analysis(f0 * (1 + 0.3e-6 * rng.standard_normal(3)))
    tab = s.sens_table(A)
    loads = []; Mb = A.I.copy(); V0 = 0.0; Tcum = 1.0
    fb, Vop = s.sweep(f0, V0, 0.0, sigma_y, rng)
    out = []
    for k in range(n_land):
        M = UB if rng.random() < 0.5 else PUB
        x = L * (0.5 + rng.uniform(-window, window)); loads.append((M * Da, x))
        V0 += sigma_V * rng.standard_normal()
        tau = sigma_tau * rng.standard_normal(); Tcum *= (1 + tau)
        if rng.random() < p_bg: loads.append((rng.uniform(18, 44) * Da, rng.uniform(0, L)))
        m = Bt['m'].copy()
        for dm, xx in loads: m += dm * s.footprint(Bt['xn'], xx)
        fa, Vop = s.sweep(c.freqs2(Bt, c.EI0, c.T0 * Tcum, m), V0, Vop, sigma_y, rng)
        r = np.log(fa / fb); fb = fa
        me, xe, te, chi = s.invert_fixed_tension(A, Mb, 0.0, r, tab)
        Mb = Mb + A.Mload(NOM * Da, xe * L)
        out.append(dict(k=k, M=M, est=float(me), xt=float(x / L), xe=float(xe), tau=float(tau), te=0.0, chi=float(chi)))
    return dict(seed=seed, truth=par, events=out)

if __name__ == '__main__':
    tag = sys.argv[1]; seeds = [int(a) for a in sys.argv[2].split(',')]
    sigma_tau = float(sys.argv[3]) if len(sys.argv) > 3 else 3e-5
    window = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    res = []; t = time.time()
    for sd in seeds:
        res.append(tube(sd, sigma_tau=sigma_tau, window=window))
        json.dump(res, open(f'random_{tag}.json', 'w'))
        print(f'{tag} {sd} {time.time() - t:.0f}s', flush=True)
