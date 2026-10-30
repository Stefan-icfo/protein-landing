"""Batch: error of central landings vs. landing index on realistic tubes (fixed-tension inversion)."""
import numpy as np, json, sys, time
import nanotube_balance_mc as s, calibrated_mc as c, final_sim as fsim

Da, L = s.Da, s.L
UB, PUB, NOM = 8564.8, 8644.8, 8604.8

def tube(seed, kind='mild', n_land=30, sigma_y=2e-6, sigma_V=0.1, p_bg=0.01, window=None):
    rng = np.random.default_rng(seed)
    Bt = fsim.make_truth(kind)
    f0 = c.freqs2(Bt, c.EI0, c.T0, Bt['m'])
    A = s.Analysis(f0 * (1 + 0.3e-6 * rng.standard_normal(3)))
    tab = s.sens_table(A)
    # per-tube scale: response to a known reference mass landing at the centre of the clean tube
    # (stands in for the in-situ calibration; its statistical error is budgeted separately)
    mref = Bt['m'] + UB * Da * s.footprint(Bt['xn'], 0.5 * L)
    rref = np.log(c.freqs2(Bt, c.EI0, c.T0, mref) / f0)
    mr, _, _, _ = s.invert_fixed_tension(A, A.I.copy(), 0.0, rref, tab)
    scale = UB / mr
    loads = []; Mb = A.I.copy(); V0 = 0.0
    fb, Vop = s.sweep(f0, V0, 0.0, sigma_y, rng)
    out = []
    for k in range(n_land):
        M = UB if rng.random() < 0.5 else PUB
        x = rng.uniform(0, L) if window is None else L * (0.5 + rng.uniform(-window, window))
        loads.append((M * Da, x))
        V0 += sigma_V * rng.standard_normal()
        if rng.random() < p_bg:
            loads.append((rng.uniform(18, 44) * Da, rng.uniform(0, L)))
        m = Bt['m'].copy()
        for dm, xx in loads:
            m += dm * s.footprint(Bt['xn'], xx)
        fa, Vop = s.sweep(c.freqs2(Bt, c.EI0, c.T0, m), V0, Vop, sigma_y, rng)
        r = np.log(fa / fb); fb = fa
        me, xe, _, chi = s.invert_fixed_tension(A, Mb, 0.0, r, tab)
        Mb = Mb + A.Mload(NOM * Da, xe * L)
        out.append(dict(k=k, M=M, est=float(me * scale), xt=float(x / L), xe=float(xe), chi=float(chi)))
    return dict(seed=seed, kind=kind, scale=float(scale), events=out)

if __name__ == '__main__':
    kind = sys.argv[1]; seeds = [int(a) for a in sys.argv[2].split(',')]; tag = sys.argv[3]
    window = float(sys.argv[4]) if len(sys.argv) > 4 else None
    res = []; t = time.time()
    for sd in seeds:
        res.append(tube(sd, kind=kind, window=window))
        json.dump(res, open(f'batch_{tag}.json', 'w'))
        print(f'{kind} {sd} {time.time() - t:.0f}s', flush=True)
