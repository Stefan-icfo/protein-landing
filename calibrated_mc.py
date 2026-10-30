"""
Calibrated version of the protocol simulation.

Truth tube imperfections (not in the analysis parametrisation):
  compliant contact: bending rigidity x0.5 over the first 3 % of the length
  2 % linear mass/EI gradient
  a 5 kDa contamination lump at 0.7 L
Per-tube calibration (an early calibration idea): 40 neutral C60 landings (720.64 Da, no charge) at random positions.
  Analysis parametrisation fitted to them: mass scale, end-segment (2 %) rigidity factors at both
  ends, linear gradient; bending rigidity and tension refitted to the unloaded frequencies each time.
Proteins are then inverted with the calibrated model (asymmetric, full-length position search).
"""
import numpy as np, sys, json, time
from scipy.linalg import eig_banded, eigh
from scipy.optimize import least_squares, minimize
import nanotube_balance_mc as b

Da, L, N = b.Da, b.L, b.N_FD
mu0, EI0, T0 = b.mu0, b.EI0, b.T0
SOFT = 0.02
C60 = 720.64

def unit_build(s_m, kL, kR, g, softL=SOFT, softR=SOFT):
    prof = lambda x: 1 + g * (2 * x / L - 1)
    def EI_fn(x):
        f = prof(x)**3
        f = np.where(x < softL * L, f * kL, f)
        f = np.where(x > (1 - softR) * L, f * kR, f)
        return f                       # unit bending rigidity profile
    return b.build(N, lambda x: mu0 * s_m * prof(x), EI_fn)

def freqs2(B, EI, T, masses, n=3, vecs=False):
    (b0, b1, b2), (t0, t1, _) = B['b'], B['t']
    s = 1 / np.sqrt(masses)
    a = np.zeros((3, len(masses)))
    a[0] = (EI * b0 + T * t0) * s * s
    a[1, :-1] = (EI * b1 + T * t1) * s[1:] * s[:-1]
    a[2, :-2] = EI * b2 * s[2:] * s[:-2]
    if vecs:
        lam, v = eig_banded(a, lower=True, select='i', select_range=(0, n - 1))
        return lam, v * s[:, None]
    lam = eig_banded(a, lower=True, select='i', select_range=(0, n - 1), eigvals_only=True)
    return np.sqrt(lam) / (2 * np.pi)

def make_truth():
    B = unit_build(1.0, 0.5, 1.0, 0.02, softL=0.03, softR=0.0)
    w = np.exp(-0.5 * ((B['xn'] - 0.7 * L) / 3e-9)**2)
    B['m'] = B['m'] + 5000 * Da * w / w.sum()
    return B

class Model:
    """Analysis model for parameters theta = (ln s_m, ln kL, ln kR, g); EI, T fitted to f0."""
    def __init__(self, theta, f0_meas, nbasis=20):
        s_m, kL, kR, g = np.exp(theta[0]), np.exp(theta[1]), np.exp(theta[2]), theta[3]
        B = unit_build(s_m, kL, kR, g)
        def res(p):
            return np.log(freqs2(B, EI0 * np.exp(p[0]), T0 * np.exp(p[1]), B['m']) / f0_meas) / 1e-6
        sol = least_squares(res, [0.0, 0.0], x_scale=[0.1, 0.1])
        self.EI, self.T = EI0 * np.exp(sol.x[0]), T0 * np.exp(sol.x[1])
        self.misfit0 = float(np.sum(sol.fun**2))
        lam, Phi = freqs2(B, self.EI, self.T, B['m'], n=nbasis, vecs=True)
        self.w2, self.Phi, self.xn = lam, Phi, B['xn']
        self.KtP = self.T * (Phi.T @ (B['Kt'] @ Phi))
        self.I = np.eye(nbasis)
        self.B = B
    def Mload(self, dm, x):
        w = b.footprint(self.xn, x); idx = w > 1e-12
        P = self.Phi[idx]
        return dm * (P.T * w[idx]) @ P
    def freqs(self, Mmod, c):
        lam = eigh(np.diag(self.w2) + c * self.KtP, Mmod, eigvals_only=True, subset_by_index=[0, 2])
        return np.sqrt(lam) / (2 * np.pi)

XG = np.linspace(0.04, 0.96, 185)

def locate(model, Mb, c_cum, r, mass_da, sig=3e-6, iters=3):
    """Known mass: find position and tension change; exact responses in the modal basis."""
    fb = model.freqs(Mb, c_cum)
    # first-order grid using the current loaded modes (from the modal problem)
    lam, U = eigh(np.diag(model.w2) + c_cum * model.KtP, Mb, subset_by_index=[0, 2])
    phi = model.Phi @ U                                  # current mass-normalised shapes at nodes
    kappa = np.array([(U[:, n] @ model.KtP @ U[:, n]) / lam[n] for n in range(3)])
    best = None
    for xf in XG:
        w = b.footprint(model.xn, xf * L)
        resp = -0.5 * mass_da * Da * (w @ phi**2)
        tau = np.dot(r - resp, kappa) / np.dot(kappa, kappa) * 2
        rr = np.sum((r - resp - 0.5 * tau * kappa)**2)
        if best is None or rr < best[0]:
            best = (rr, xf, tau)
    _, x, tau = best
    p = np.array([x, tau])
    def g(q):
        fa = model.freqs(Mb + model.Mload(mass_da * Da, q[0] * L), c_cum + q[1])
        return np.log(fa / fb)
    for _ in range(iters):                               # Gauss-Newton on (x, tau), exact model
        g0 = g(p)
        J = np.column_stack([(g(p + [1e-4, 0]) - g0) / 1e-4, (g(p + [0, 1e-6]) - g0) / 1e-6])
        dp, *_ = np.linalg.lstsq(J, r - g0, rcond=None)
        p = p + dp; p[0] = np.clip(p[0], 0.02, 0.98)
    res = (g(p) - r) / sig
    return p, float(res @ res)

def calib_cost(theta, f0_meas, cal_r):
    try:
        model = Model(theta, f0_meas)
    except Exception:
        return 1e12
    Mb = model.I.copy(); c = 0.0; tot = 0.0
    for r in cal_r:
        (x, tau), chi = locate(model, Mb, c, r, C60)
        tot += chi
        Mb = Mb + model.Mload(C60 * Da, x * L); c += tau
    return tot

def invert_unknown(model, Mb, c_cum, r, sig=3e-6):
    """Unknown mass: scan position (full length), inner Gauss-Newton on (mass, tension)."""
    fb = model.freqs(Mb, c_cum)
    def g(m_da, xf, tau):
        return np.log(model.freqs(Mb + model.Mload(m_da * Da, xf * L), c_cum + tau) / fb)
    def inner(xf, m0, t0, iters):
        p = np.array([m0, t0])
        try:
            for _ in range(iters):
                g0 = g(p[0], xf, p[1])
                J = np.column_stack([g(p[0] + 1.0, xf, p[1]) - g0, (g(p[0], xf, p[1] + 1e-6) - g0) / 1e-6])
                dp, *_ = np.linalg.lstsq(J, r - g0, rcond=None)
                p = p + dp; p[0] = max(p[0], 0.0)
            rr = g(p[0], xf, p[1]) - r
            return p, float(rr @ rr)
        except np.linalg.LinAlgError:
            return p, np.inf
    # first-order start per grid point from unloaded-state sensitivities
    scan = []
    for xf in XG:
        p, chi = inner(xf, 8000.0, 0.0, 2)
        scan.append((chi, xf, p))
    chi, xb, pb = min(scan, key=lambda t: t[0])
    fine = []
    for xf in np.clip(xb + np.linspace(-0.006, 0.006, 25), 0.04, 0.96):
        p, chi = inner(xf, pb[0], pb[1], 3)
        fine.append((chi, xf, p))
    j = int(np.argmin([t[0] for t in fine]))
    if 0 < j < len(fine) - 1:
        (c0, x0, _), (c1, x1, _), (c2, x2, _) = fine[j - 1], fine[j], fine[j + 1]
        den = c0 - 2 * c1 + c2
        xh = x1 + 0.5 * (x0 - x2) * (c0 - c2) / (2 * den) if den > 0 else x1
        xh = float(np.clip(xh, min(x0, x2), max(x0, x2)))
    else:
        xh = fine[j][1]
    p, chi = inner(xh, fine[j][2][0], fine[j][2][1], 3)
    return p[0], xh, p[1], chi / sig**2

def run_tube(seed, cap=0.3, n_cal=40, sigma_y=2e-6, sigma_tau=3e-4, sigma_V=0.1, p_bg=0.01,
             proteins=(8564.8, 8644.8), calibrate=True):
    rng = np.random.default_rng(seed)
    Bt = make_truth()
    m_true_tube = Bt['m'].sum()
    f0 = freqs2(Bt, EI0, T0, Bt['m'])
    f0_meas = f0 * (1 + 0.3e-6 * rng.standard_normal(3))
    true_loads = []; Tcum = 1.0; V0 = 0.0
    fb, Vop = b.sweep(f0, V0, 0.0, sigma_y, rng)
    def land(mass_da, charged):
        nonlocal Tcum, V0, fb, Vop
        x = rng.uniform(0, L)
        tau = sigma_tau * rng.standard_normal(); Tcum *= (1 + tau)
        if charged: V0 += sigma_V * rng.standard_normal()
        true_loads.append((mass_da * Da, x))
        bg = 0.0
        if rng.random() < p_bg:
            bg = rng.uniform(18, 44); true_loads.append((bg * Da, rng.uniform(0, L)))
        m = Bt['m'].copy()
        for dm, xx in true_loads:
            m += dm * b.footprint(Bt['xn'], xx)
        fm = freqs2(Bt, EI0, T0 * Tcum, m)
        fa, Vop_new = b.sweep(fm, V0, Vop, sigma_y, rng)
        r = np.log(fa / fb)
        fb, Vop = fa, Vop_new
        return r, x / L, bg
    # --- calibration series with neutral C60
    cal_r = [land(C60, False)[0] for _ in range(n_cal)]
    theta0 = np.zeros(4)
    if calibrate:
        opt = minimize(calib_cost, theta0, args=(f0_meas, cal_r), method='Nelder-Mead',
                       options=dict(xatol=1e-4, fatol=0.5, maxiter=400, initial_simplex=np.array(
                           [[0, 0, 0, 0], [0.02, 0, 0, 0], [0, -0.3, 0, 0], [0, 0, -0.3, 0], [0, 0, 0, 0.02]])))
        theta = opt.x
    else:
        theta = theta0
    model = Model(theta, f0_meas)
    # replay calibration landings to place them in the model
    Mb = model.I.copy(); c = 0.0
    for r in cal_r:
        (x, tau), _ = locate(model, Mb, c, r, C60)
        Mb = Mb + model.Mload(C60 * Da, x * L); c += tau
    # --- protein landings
    out = []
    added = 0.0
    while added < cap * m_true_tube:
        M = proteins[rng.integers(len(proteins))]
        load_before = added / m_true_tube
        r, xt, bg = land(M, True)
        added += M * Da
        m_e, x_e, t_e, chi = invert_unknown(model, Mb, c, r)
        Mb = Mb + model.Mload(m_e * Da, x_e * L); c += t_e
        out.append(dict(M=M, est=float(m_e), xt=float(xt), xe=float(x_e), load=load_before,
                        chi=chi, bg=bg, theta=[float(v) for v in theta]))
    return out

if __name__ == '__main__':
    seeds = [int(a) for a in sys.argv[1].split(',')]
    cal = sys.argv[2] == 'cal'
    tag = sys.argv[3]
    res = []
    t = time.time()
    for sd in seeds:
        res += [dict(e, tube=sd) for e in run_tube(sd, calibrate=cal)]
        print(f'seed {sd} done {time.time() - t:.0f}s', flush=True)
    json.dump(res, open(f'out_{tag}.json', 'w'))
