"""
Protocol Monte Carlo with C60 scale calibration (earlier version; make_truth is still used).

Truth tube (unknown to the analysis):
  'mild'  : compliant contact (bending rigidity x0.8 over 3 % of the length) at one end,
            1 % linear mass/rigidity gradient, 2 kDa contamination lump at 0.7 L
  'ideal' : uniform, symmetric clamped beam (mode shapes known)
Analysis: uniform symmetric clamped beam, EI and T fitted to the unloaded frequencies;
  three modes; five-point gate sweeps; inversion for (mass, position, tension change).
Per tube:
  1) calibration: n_cal neutral C60 landings (720.64 Da) -> per-tube scale factor from the landings
     estimated within 0.1 L of the centre (no protein involved);
  2) proteins: ubiquitin (8564.8 Da) or Ser65-phospho-ubiquitin (8644.8 Da), 50:50, until the
     added load reaches `cap` x tube mass. Non-central landings are book-kept at the nominal
     selected mass and their estimated position.
Nuisances per landing: tension jump N(0, sigma_tau), gate-vertex (charge) shift N(0, 0.1 V) for
  proteins, 1 % background molecules (18-44 Da) in the same step, frequency noise sigma_y.
"""
import numpy as np, json, sys, time
import nanotube_balance_mc as s
import calibrated_mc as c

Da, L = s.Da, s.L
C60 = 720.64
UB, PUB = 8564.8, 8644.8

def make_truth(kind):
    if kind == 'ideal':
        B = c.unit_build(1.0, 1.0, 1.0, 0.0)
    else:
        B = c.unit_build(1.0, 0.8, 1.0, 0.01, softL=0.03, softR=0.0)
        w = np.exp(-0.5 * ((B['xn'] - 0.7 * L) / 3e-9)**2)
        B['m'] = B['m'] + 2000 * Da * w / w.sum()
    return B

def run_tube(seed, kind='mild', n_cal=120, cap=0.3, sigma_y=2e-6, sigma_tau=3e-4, sigma_V=0.1, p_bg=0.01):
    rng = np.random.default_rng(seed)
    Bt = make_truth(kind)
    m_tube = Bt['m'].sum()
    f0 = c.freqs2(Bt, c.EI0, c.T0, Bt['m'])
    A = s.Analysis(f0 * (1 + 0.3e-6 * rng.standard_normal(3)))
    tab = s.sens_table(A)
    Mb = A.I.copy(); cc = 0.0
    loads = []; Tcum = 1.0; V0 = 0.0
    state = dict(fb=None, Vop=0.0)
    state['fb'], state['Vop'] = s.sweep(f0, V0, 0.0, sigma_y, rng)
    def land(mass_da, charged):
        nonlocal Tcum, V0
        x = rng.uniform(0, L)
        Tcum *= (1 + sigma_tau * rng.standard_normal())
        if charged:
            V0 += sigma_V * rng.standard_normal()
        loads.append((mass_da * Da, x))
        if rng.random() < p_bg:
            loads.append((rng.uniform(18, 44) * Da, rng.uniform(0, L)))
        m = Bt['m'].copy()
        for dm, xx in loads:
            m += dm * s.footprint(Bt['xn'], xx)
        fa, Vn = s.sweep(c.freqs2(Bt, c.EI0, c.T0 * Tcum, m), V0, state['Vop'], sigma_y, rng)
        r = np.log(fa / state['fb'])
        state['fb'], state['Vop'] = fa, Vn
        xf = x / L
        return r, (xf if xf <= 0.5 else 1 - xf)
    # 1) calibration with neutral C60
    ratios = []
    for _ in range(n_cal):
        r, xt = land(C60, False)
        m_e, x_e, t_e, chi = s.invert(A, Mb, cc, r, tab)
        if abs(x_e - 0.5) <= 0.1 and m_e > 0:
            ratios.append(C60 / m_e)
        Mb = Mb + A.Mload(C60 * Da, x_e * L); cc += t_e
    ratios = np.array(ratios)
    med = np.median(ratios); mad = 1.4826 * np.median(np.abs(ratios - med))
    keep = np.abs(ratios - med) <= 3 * mad if mad > 0 else np.ones_like(ratios, bool)
    scale = float(ratios[keep].mean())
    scale_err = float(ratios[keep].std(ddof=1) / np.sqrt(keep.sum()))
    # 2) proteins
    out = []
    added = sum(dm for dm, _ in loads)
    while added < cap * m_tube:
        M = UB if rng.random() < 0.5 else PUB
        load_before = added / m_tube
        r, xt = land(M, True)
        added += M * Da
        m_e, x_e, t_e, chi = s.invert(A, Mb, cc, r, tab)
        Mb = Mb + A.Mload(8604.8 * Da, x_e * L); cc += t_e
        out.append(dict(M=M, raw=float(m_e), est=float(m_e * scale), xt=float(xt), xe=float(x_e),
                        load=float(load_before), chi=float(chi)))
    return dict(seed=seed, kind=kind, scale=scale, scale_err=scale_err, n_cal_acc=int(keep.sum()), events=out)

if __name__ == '__main__':
    kind = sys.argv[1]
    seeds = [int(a) for a in sys.argv[2].split(',')]
    tag = sys.argv[3]
    res = []; t = time.time()
    for sd in seeds:
        res.append(run_tube(sd, kind=kind))
        print(f'{kind} seed {sd} done {time.time() - t:.0f}s  scale {res[-1]["scale"]:.5f} +- {res[-1]["scale_err"]:.5f}', flush=True)
        json.dump(res, open(f'final_{tag}.json', 'w'))
