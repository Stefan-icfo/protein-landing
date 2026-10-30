"""
Monte Carlo of the event-resolved nanotube weighing protocol.

Truth ("the real tube"): clamped-clamped beam under tension, finite differences (N=399),
  with two features the analysis does NOT know about:
  - 3 % linear diameter gradient (mass per length and bending rigidity vary along the tube)
  - a compliant contact: bending rigidity halved over the first 2 % of the length.
Analysis model: uniform, symmetric clamped beam; tube mass assumed known (scale calibrated
  separately); bending rigidity and tension fitted to the measured unloaded frequencies of the
  three modes; loaded responses computed in a 20-mode modal basis (well conditioned, smooth).
Landed protein: Gaussian mass footprint, 1.5 nm (1 sigma), uniformly random position.
Per landing nuisances:
  - tension change      tau ~ N(0, sigma_tau)  (relative)
  - charge offset       vertex of f(V_gate) moves by N(0, sigma_V)
  - background molecule with probability p_bg (18-44 Da, random position, same step)
  - frequency noise     sigma_y per reading (Allan deviation at 1 s)
Readout: three flexural modes; before and after each landing, a 5-point gate sweep (+-0.4 V around
  the last estimated vertex); vertex frequency from a parabola fit (truth has a quartic term too).
Inversion: (mass, position, tension change) from the three mode-frequency ratios, with all previously
  estimated loads included in the analysis model.
Landings continue on a tube until the true added load reaches `cap` x tube mass.
"""
import numpy as np, time
from scipy.linalg import eig_banded, eigh
from scipy.optimize import least_squares

Da = 1.66053906660e-27

# ---------------- tube ----------------
L = 300e-9
d0 = 2e-9
mu0 = 7.6e-7 * np.pi * d0                  # graphene areal density x circumference, kg/m
EI0 = 1e12 * np.pi * d0**3 * 0.34e-9 / 8   # E = 1 TPa, wall 0.34 nm
T0 = 0.5e-9                                # N
m_tube = mu0 * L
N_FD = 399
FOOT = 1.5e-9                              # 1-sigma footprint of a landed protein

def build(N, mu_fn, EI_fn, clamp=(1.0, 1.0)):
    h = L / (N + 1)
    xn = np.arange(N + 2) * h
    mu = mu_fn(xn); EI = EI_fn(xn)
    C = np.zeros((N + 2, N))
    for j in range(N + 2):
        for k, c in ((j - 1, 1.0), (j, -2.0), (j + 1, 1.0)):
            if k == -1:
                k, c = 1, c * clamp[0]
            if k == N + 2:
                k, c = N, c * clamp[1]
            if 1 <= k <= N:
                C[j, k - 1] += c / h**2
    w = np.ones(N + 2); w[0] = w[-1] = 0.5
    Kb = h * C.T @ (np.diag(EI * w) @ C)
    S = np.zeros((N + 1, N))
    for j in range(N + 1):
        if 1 <= j + 1 <= N: S[j, j] += 1 / h
        if 1 <= j <= N: S[j, j - 1] -= 1 / h
    Kt = h * S.T @ S
    bands = lambda K: (np.diag(K).copy(), np.diag(K, -1).copy(), np.diag(K, -2).copy())
    return dict(N=N, h=h, xn=xn[1:N + 1], b=bands(Kb), t=bands(Kt), Kt=Kt, m=mu[1:N + 1] * h)

def _banded_A(B, T, masses):
    (b0, b1, b2), (t0, t1, _) = B['b'], B['t']
    s = 1 / np.sqrt(masses)
    a = np.zeros((3, len(masses)))
    a[0] = (b0 + T * t0) * s * s
    a[1, :-1] = (b1 + T * t1) * s[1:] * s[:-1]
    a[2, :-2] = b2 * s[2:] * s[:-2]
    return a, s

def freqs(B, T, masses, n=3):
    a, _ = _banded_A(B, T, masses)
    lam = eig_banded(a, lower=True, select='i', select_range=(0, n - 1), eigvals_only=True)
    return np.sqrt(lam) / (2 * np.pi)

def footprint(xn, x):
    w = np.exp(-0.5 * ((xn - x) / FOOT)**2)
    return w / w.sum() if w.sum() > 0 else w

def loaded(B, loads):
    m = B['m'].copy()
    for dm, x in loads:
        m += dm * footprint(B['xn'], x)
    return m

def make_truth(grad=0.03, soft_len=0.02, soft_fac=0.5):
    prof = lambda x: 1 + grad * (2 * x / L - 1)
    return build(N_FD, lambda x: mu0 * prof(x),
                 lambda x: EI0 * prof(x)**3 * np.where(x < soft_len * L, soft_fac, 1.0))

class Analysis:
    """Uniform symmetric clamped beam, tube mass known, EI and T fitted; modal basis for loads."""
    def __init__(self, f_meas, nbasis=20):
        def model(p):
            EI, T = EI0 * np.exp(p[0]), T0 * np.exp(p[1])
            return build(N_FD, lambda x: mu0 * np.ones_like(x), lambda x: EI * np.ones_like(x)), T
        def res(p):
            B, T = model(p)
            return np.log(freqs(B, T, B['m']) / f_meas) / 1e-6
        sol = least_squares(res, [0.0, 0.0], x_scale=[0.1, 0.1])
        B, T = model(sol.x)
        a, s = _banded_A(B, T, B['m'])
        lam, v = eig_banded(a, lower=True, select='i', select_range=(0, nbasis - 1))
        self.Phi = v * s[:, None]                 # mass-normalised mode shapes
        self.w2 = lam
        self.KtP = T * (self.Phi.T @ (B['Kt'] @ self.Phi))
        self.xn = B['xn']; self.K = nbasis
        self.I = np.eye(nbasis)
    def Mload(self, dm, x):
        w = footprint(self.xn, x)
        idx = w > 1e-12
        P = self.Phi[idx]
        return dm * (P.T * w[idx]) @ P
    def freqs(self, Mmod, c):
        lam = eigh(np.diag(self.w2) + c * self.KtP, Mmod, eigvals_only=True, subset_by_index=[0, 2])
        return np.sqrt(lam) / (2 * np.pi)

XS = np.linspace(0.04, 0.96, 185)

def sens_table(A, probe_da=1000.0):
    f0 = A.freqs(A.I, 0.0)
    S = np.array([np.log(A.freqs(A.I + A.Mload(probe_da * Da, xf * L), 0.0) / f0) / probe_da for xf in XS])
    st = np.log(A.freqs(A.I, 1e-4) / f0) / 1e-4
    return S, st

def _inner(A, Mb, c_cum, fb, r, xf, dm0, tau0, iters=3):
    """Gauss-Newton for (mass in Da, tension change) at fixed position; returns params and chi2."""
    p = np.array([dm0, tau0], float)
    def model(q):
        fa = A.freqs(Mb + A.Mload(q[0] * Da, xf * L), c_cum + q[1])
        return np.log(fa / fb)
    p[0] = max(p[0], 0.0)
    try:
        for _ in range(iters):
            g0 = model(p)
            J = np.column_stack([(model(p + [1.0, 0]) - g0) / 1.0, (model(p + [0, 1e-6]) - g0) / 1e-6])
            dp, *_ = np.linalg.lstsq(J, r - g0, rcond=None)
            p = p + dp
            p[0] = max(p[0], 0.0)
        resid = model(p) - r
        return p, float(resid @ resid)
    except np.linalg.LinAlgError:
        return p, np.inf

def invert(A, Mb, c_cum, r_meas, table, sig=3e-6):
    S, st = table
    fb = A.freqs(Mb, c_cum)
    scan = []
    for i, xf in enumerate(XS):                      # coarse scan, nonlinear inner solve
        G = np.column_stack([S[i], st])
        p0, *_ = np.linalg.lstsq(G, r_meas, rcond=None)
        p, chi = _inner(A, Mb, c_cum, fb, r_meas, xf, p0[0], p0[1], iters=2)
        scan.append((chi, xf, p))
    chi, xb, pb = min(scan, key=lambda t: t[0])
    fine = []
    for xf in np.clip(xb + np.linspace(-0.006, 0.006, 25), 0.04, 0.96):
        p, chi = _inner(A, Mb, c_cum, fb, r_meas, xf, pb[0], pb[1], iters=3)
        fine.append((chi, xf, p))
    fine.sort(key=lambda t: t[1])
    j = int(np.argmin([t[0] for t in fine]))
    if 0 < j < len(fine) - 1:                        # parabolic refinement of the position
        (c0, x0, _), (c1, x1, _), (c2, x2, _) = fine[j - 1], fine[j], fine[j + 1]
        den = (c0 - 2 * c1 + c2)
        xh = x1 + 0.5 * (x0 - x2) * (c0 - c2) / (2 * den) if den > 0 else x1
        xh = float(np.clip(xh, min(x0, x2), max(x0, x2)))
    else:
        xh = fine[j][1]
    p, chi = _inner(A, Mb, c_cum, fb, r_meas, xh, fine[j][2][0], fine[j][2][1], iters=3)
    return p[0], xh, p[1], chi / sig**2

# gate response (relative frequency vs gate voltage around the flat point)
A_GATE = np.array([2e-3, 1.5e-3, 1.2e-3])    # 1/V^2
B_GATE = np.array([4e-5, 3e-5, 2e-5])        # 1/V^4 (not in the analysis fit)
DV = np.array([-0.4, -0.2, 0.0, 0.2, 0.4])

def sweep(fmech, V0_true, Vop, sigma_y, rng):
    V = Vop + DV
    est_f = np.empty(3); vertices = []
    for n in range(3):
        u = V - V0_true
        f = fmech[n] * (1 + A_GATE[n] * u**2 + B_GATE[n] * u**4) * (1 + sigma_y * rng.standard_normal(len(V)))
        c2, c1, c0 = np.polyfit(V - Vop, f / fmech[n], 2)
        est_f[n] = fmech[n] * (c0 - c1**2 / (4 * c2))
        vertices.append(Vop - c1 / (2 * c2))
    return est_f, float(np.mean(vertices))

def run(n_tubes, cap, sigma_y=2e-6, sigma_tau=3e-4, sigma_V=0.1, p_bg=0.01,
        mismatch=True, masses_da=(8564.8,), seed=0, verbose=False):
    rng = np.random.default_rng(seed)
    Bt = make_truth() if mismatch else make_truth(grad=0.0, soft_len=0.0)
    f0 = freqs(Bt, T0, Bt['m'])
    out = []
    for tube in range(n_tubes):
        true_loads = []
        Tcum_t, c_e = 1.0, 0.0
        V0_t = 0.0
        A = Analysis(f0 * (1 + 0.3e-6 * rng.standard_normal(3)))   # averaged unloaded characterisation
        table = sens_table(A)
        Mb = A.I.copy()
        fb, Vop = sweep(f0, V0_t, 0.0, sigma_y, rng)
        added = 0.0; k = 0
        while added < cap * m_tube:
            M = masses_da[rng.integers(len(masses_da))]
            dm = M * Da; x = rng.uniform(0, L)
            tau = sigma_tau * rng.standard_normal()
            V0_t += sigma_V * rng.standard_normal()
            Tcum_t *= (1 + tau)
            load_before = added / m_tube
            true_loads.append((dm, x)); added += dm
            bg = 0.0
            if rng.random() < p_bg:
                bg = rng.uniform(18, 44) * Da
                true_loads.append((bg, rng.uniform(0, L))); added += bg
            fmech = freqs(Bt, T0 * Tcum_t, loaded(Bt, true_loads))
            fa, Vop_new = sweep(fmech, V0_t, Vop, sigma_y, rng)
            dm_e, x_e, tau_e, cost = invert(A, Mb, c_e, np.log(fa / fb), table)
            Mb = Mb + A.Mload(dm_e * Da, x_e * L); c_e += tau_e
            xf = x / L if x / L <= 0.5 else 1 - x / L
            out.append(dict(tube=tube, k=k, M=M, est=dm_e, xt=xf, xe=x_e, tau=tau, tau_e=tau_e,
                            bg=bg / Da, load=load_before, cost=cost))
            fb, Vop = fa, Vop_new
            k += 1
        if verbose:
            print(f'tube {tube}: {k} landings', flush=True)
    return out

if __name__ == '__main__':
    t = time.time()
    test = run(1, 0.08, seed=1)
    print(len(test), 'events', round(time.time() - t, 1), 's')
    for e in test:
        print(round(e['xt'], 3), round(e['xe'], 3), round(e['est'] - e['M'], 1), round(e['cost'], 2))


def invert_fixed_tension(A, Mb, c_cum, r_meas, table, sig=3e-6):
    """Two unknowns (mass, position); tension held at its current estimate. Three modes give one
    redundant equation, which selects the landing side on an asymmetric (loaded) tube."""
    S, st = table
    fb = A.freqs(Mb, c_cum)
    def g(m_da, xf):
        return np.log(A.freqs(Mb + A.Mload(m_da * Da, xf * L), c_cum) / fb)
    def inner(xf, m0, iters):
        m = max(m0, 1.0)
        try:
            for _ in range(iters):
                g0 = g(m, xf); J = g(m + 1.0, xf) - g0
                m = max(m + float(J @ (r_meas - g0) / (J @ J)), 0.0)
            rr = g(m, xf) - r_meas
            return m, float(rr @ rr)
        except np.linalg.LinAlgError:
            return m, np.inf
    scan = []
    for i, xf in enumerate(XS):
        m0 = float(S[i] @ r_meas / (S[i] @ S[i]))
        m, chi = inner(xf, m0, 2)
        scan.append((chi, xf, m))
    chi, xb, mb = min(scan, key=lambda t: t[0])
    fine = []
    for xf in np.clip(xb + np.linspace(-0.006, 0.006, 25), 0.04, 0.96):
        m, chi = inner(xf, mb, 3); fine.append((chi, xf, m))
    j = int(np.argmin([t[0] for t in fine]))
    if 0 < j < len(fine) - 1:
        (c0, x0, _), (c1, x1, _), (c2, x2, _) = fine[j - 1], fine[j], fine[j + 1]
        den = c0 - 2 * c1 + c2
        xh = x1 + 0.5 * (x0 - x2) * (c0 - c2) / (2 * den) if den > 0 else x1
        xh = float(np.clip(xh, min(x0, x2), max(x0, x2)))
    else:
        xh = fine[j][1]
    m, chi = inner(xh, fine[j][2], 3)
    return m, xh, 0.0, chi / sig**2
