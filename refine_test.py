"""Ideal tube, no noise, true bookkeeping: grid-based inversion vs. a continuous least-squares refinement of (mass, position)."""
import sys, numpy as np
from scipy.optimize import least_squares
import nanotube_balance_mc as s, calibrated_mc as c, final_sim as fsim
from analyze_tension import stats
Da, L = s.Da, s.L
UB, PUB = 8564.8, 8644.8
def refine(A, Mb, r, m0, x0):
    fb = A.freqs(Mb, 0.0)
    def res(p):
        return (np.log(A.freqs(Mb + A.Mload(p[0] * Da, p[1] * L), 0.0) / fb) - r) / 1e-7
    sol = least_squares(res, [m0, x0], x_scale=[10.0, 1e-3], diff_step=[1e-5, 1e-6], xtol=1e-12, ftol=1e-12)
    return sol.x
def run(kind, refine_on):
    R = []
    for sd in range(400, 404):
        rng = np.random.default_rng(sd)
        Bt = fsim.make_truth(kind); f0 = c.freqs2(Bt, c.EI0, c.T0, Bt['m'])
        A = s.Analysis(f0); tab = s.sens_table(A)
        loads = []; Mb = A.I.copy(); fb = f0; out = []
        for k in range(30):
            M = UB if rng.random() < 0.5 else PUB
            x = L * (0.5 + rng.uniform(-0.05, 0.05)); loads.append((M * Da, x))
            m = Bt['m'].copy()
            for dm, xx in loads: m += dm * s.footprint(Bt['xn'], xx)
            fa = c.freqs2(Bt, c.EI0, c.T0, m); r = np.log(fa / fb); fb = fa
            me, xe, _, _ = s.invert_fixed_tension(A, Mb, 0.0, r, tab)
            if refine_on: me, xe = refine(A, Mb, r, me, xe)
            Mb = Mb + A.Mload(M * Da, x)
            out.append(dict(k=k, M=M, est=float(me), xe=float(xe), xt=float(x / L), chi=0.0))
        R.append(dict(seed=sd, events=out))
    o, e, M = stats(R)
    return o
for kind in sys.argv[1:]:
    for ref in (False, True):
        o = run(kind, ref)
        print(f'{kind:5s} refine={ref}: raw SD {o["raw_sd"]:.1f} Da, after per-tube correction SD {o["sd"]:.1f} Da', flush=True)
