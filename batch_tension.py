"""Slit landings (+-0.05 L) with tension changes injected at every landing.
Inversion with the tension change as a third unknown (all three modes used, no redundancy),
or with tension held fixed (for comparison). Otherwise identical to batch.py."""
import numpy as np, json, sys, time
import nanotube_balance_mc as s, calibrated_mc as c, final_sim as fsim

Da, L = s.Da, s.L
UB, PUB, NOM = 8564.8, 8644.8, 8604.8

def tube(seed, sigma_tau, free=True, n_land=30, sigma_y=2e-6, sigma_V=0.1, p_bg=0.01, window=0.05):
    rng = np.random.default_rng(seed)
    Bt = fsim.make_truth('mild')
    f0 = c.freqs2(Bt, c.EI0, c.T0, Bt['m'])
    A = s.Analysis(f0 * (1 + 0.3e-6 * rng.standard_normal(3)))
    tab = s.sens_table(A)
    loads = []; Mb = A.I.copy(); V0 = 0.0; Tcum = 1.0; cc = 0.0
    fb, Vop = s.sweep(f0, V0, 0.0, sigma_y, rng)
    out = []
    for k in range(n_land):
        M = UB if rng.random() < 0.5 else PUB
        x = L * (0.5 + rng.uniform(-window, window))
        loads.append((M * Da, x))
        V0 += sigma_V * rng.standard_normal()
        tau = sigma_tau * rng.standard_normal(); Tcum *= (1 + tau)
        if rng.random() < p_bg:
            loads.append((rng.uniform(18, 44) * Da, rng.uniform(0, L)))
        m = Bt['m'].copy()
        for dm, xx in loads:
            m += dm * s.footprint(Bt['xn'], xx)
        fa, Vop = s.sweep(c.freqs2(Bt, c.EI0, c.T0 * Tcum, m), V0, Vop, sigma_y, rng)
        r = np.log(fa / fb); fb = fa
        if free:
            me, xe, te, chi = s.invert(A, Mb, cc, r, tab)
        else:
            me, xe, te, chi = s.invert_fixed_tension(A, Mb, cc, r, tab)
        Mb = Mb + A.Mload(NOM * Da, xe * L); cc += te
        out.append(dict(k=k, M=M, est=float(me), xt=float(x / L), xe=float(xe), tau=float(tau), te=float(te), chi=float(chi)))
    return dict(seed=seed, sigma_tau=sigma_tau, free=free, events=out)

if __name__ == '__main__':
    sigma_tau = float(sys.argv[1]); free = sys.argv[2] == 'free'; tag = sys.argv[3]
    seeds = [int(a) for a in sys.argv[4].split(',')] if len(sys.argv) > 4 else range(300, 320)
    res = []; t = time.time()
    for sd in seeds:
        res.append(tube(sd, sigma_tau, free))
        json.dump(res, open(f'tension_{tag}.json', 'w'))
        print(f'{tag} {sd} {time.time() - t:.0f}s', flush=True)
