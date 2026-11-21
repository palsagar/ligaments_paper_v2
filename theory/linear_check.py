"""Linear Rayleigh-Plateau check of ligament breakup statistics.

Units: mean radius R = 1, capillary time sqrt(rho R^3 / sigma) = 1.
Usage: python linear_check.py [N_ligaments] [N_convergence_subset]
"""
import sys
import time
import numpy as np
from scipy import special, optimize

L, LAM, W, OH = 100.0, 9.02, 2.0, 0.01
SEED = 12345
P = print


# ----------------------------------------------------------------- dispersion
def weber_roots(x, oh=OH):
    x = np.asarray(x, float)
    b = 3 * oh * x**2
    c = -0.5 * x**2 * (1 - x**2)
    d = np.sqrt((b * b - 4 * c).astype(complex))
    return (-b + d) / 2, (-b - d) / 2


def growth(x, t):
    """G[t, mode]: amplitude factor of a mode that starts with zero velocity."""
    s1, s2 = weber_roots(x)
    t = np.asarray(t, float)[:, None]
    return ((s1 * np.exp(s2 * t) - s2 * np.exp(s1 * t)) / (s1 - s2)).real


def exact_rate(x, oh):
    """Growth rate from the exact viscous dispersion relation of a liquid
    cylinder (Rayleigh 1892, Chandrasekhar 1961), with y^2 = x^2 + s/Oh."""
    I0x, I1x, dI1x = special.iv(0, x), special.iv(1, x), special.ivp(1, x)

    def f(s):
        y = np.sqrt(x * x + s / oh)
        ry = special.ive(0, y) / special.ive(1, y) - 1 / y  # I1'(y)/I1(y)
        lhs = (2 * x * x * (x * x + y * y) * (dI1x / I0x)
               * (1 - (2 * x * y / (x * x + y * y)) * (I1x / dI1x) * ry)
               - (x**4 - y**4))
        return lhs - x * (1 - x * x) * I1x / I0x / oh**2

    return optimize.brentq(f, 1e-9, 5.0, xtol=1e-13)


def dispersion_checks():
    P("=" * 78)
    P("DISPERSION RELATION CHECKS")
    x = np.linspace(0.001, 0.999, 999)
    sw = weber_roots(x)[0].real
    i = sw.argmax()
    P(f"Weber, Oh=0.01: max growth rate {sw[i]:.4f} at x = {x[i]:.3f} "
      f"(wavelength {2*np.pi/x[i]:.3f})")
    # validation of the exact relation in two known limits
    si = np.array([exact_rate(xx, 1e-6) for xx in x])
    j = si.argmax()
    ray = np.sqrt(x * (1 - x**2) * special.iv(1, x) / special.iv(0, x))
    P(f"Exact relation, Oh=1e-6: max {si[j]:.4f} at x = {x[j]:.3f}; "
      f"max |s - Rayleigh inviscid| = {np.abs(si-ray).max():.2e}  (expected 0.3433 at 0.697)")
    oh, xs = 30.0, x[x > 0.3]
    sv = np.array([exact_rate(xx, oh) for xx in xs])
    r = special.iv(0, xs) / special.iv(1, xs)
    stokes = (1 - xs**2) / (2 * oh * (xs**2 * r**2 - 1 - xs**2))
    P(f"Exact relation, Oh=30, 0.3<x<1: max relative difference from Rayleigh's Stokes-limit "
      f"formula = {np.abs(sv/stokes-1).max():.2e}; from Weber = "
      f"{np.abs(sv/weber_roots(xs, oh)[0].real-1).max():.2e}")
    se = np.array([exact_rate(xx, OH) for xx in x])
    j = se.argmax()
    P(f"Exact relation, Oh=0.01: max {se[j]:.4f} at x = {x[j]:.3f} "
      f"(wavelength {2*np.pi/x[j]:.3f})")
    d = sw - se
    P(f"Weber minus exact over 0<x<1: max abs {np.abs(d).max():.4f} at x={x[np.abs(d).argmax()]:.3f}; "
      f"at the peak {sw.max()-se.max():+.4f} ({(sw.max()/se.max()-1)*100:+.2f}%)")
    for xx in (0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95):
        k = np.abs(x - xx).argmin()
        P(f"   x={xx:4.2f}: Weber {sw[k]:.4f}  exact {se[k]:.4f}  diff {d[k]:+.4f} "
          f"({(sw[k]/se[k]-1)*100:+.1f}%)")
    n = np.arange(1, 16)
    xn = 2 * np.pi * n / L
    P("Box modes n=1..15 (x_n, Weber rate, exact rate):")
    P("   " + "  ".join(f"{a}:{b:.3f}/{c:.4f}/{exact_rate(b, OH):.4f}"
                        for a, b, c in zip(n, xn, weber_roots(xn)[0].real)))
    # curvature of the Weber growth curve at its maximum (for the bandwidth estimate)
    h = 1e-3
    xm = x[i]
    s2 = (weber_roots(xm + h)[0].real - 2 * weber_roots(xm)[0].real
          + weber_roots(xm - h)[0].real) / h**2
    return xm, sw[i], s2


# ----------------------------------------------------------------- pinch sites
def find_sites(A, B, nx=2048, dt=0.02, r=10, tmax=40.0):
    """A, B: (N, nmodes) cosine and sine coefficients of modes n = 1..nmodes.
    tau(x) = first time with 1 + h <= 0. It is bracketed on a coarse time grid
    (step r*dt), then located on the fine grid (step dt) and linearly
    interpolated. r = 1 is brute force on the fine grid.
    Returns flat arrays (ligament index, x_i, tau_i), sorted by ligament and x."""
    N, nm = A.shape
    k = 2 * np.pi * np.arange(1, nm + 1) / L
    dx = L / nx
    x = np.arange(nx) * dx
    COS, SIN = np.cos(np.outer(k, x)), np.sin(np.outer(k, x))
    ntc = int(round(tmax / (dt * r)))
    tf = np.arange(ntc * r + 1) * dt
    Gf = growth(k, tf)
    Gc = Gf[::r]
    bs = max(1, min(int(2.5e7 / ((ntc + 1) * nx)), int(2.5e7 / (nx * (r + 1) * nm))))
    out = []
    for s in range(0, N, bs):
        C = A[s:s + bs, :, None] * COS + B[s:s + bs, :, None] * SIN  # (b, nm, nx)
        m = (Gc @ C) <= -1.0                                         # (b, ntc+1, nx)
        j = m.argmax(axis=1)
        hit = m.any(axis=1)
        del m
        tau = np.full(j.shape, np.inf)
        tau[hit & (j == 0)] = 0.0
        bi, xi = np.nonzero(hit & (j > 0))
        fi = ((j[bi, xi] - 1) * r)[:, None] + np.arange(r + 1)
        hf = np.einsum('mtn,mn->mt', Gf[fi], C[bi, :, xi])
        mf = hf <= -1.0
        q = np.where(mf.any(axis=1), mf.argmax(axis=1), r).clip(1)
        ar = np.arange(q.size)
        h0, h1 = hf[ar, q - 1], hf[ar, q]
        tau[bi, xi] = tf[fi[ar, q - 1]] + dt * (h0 + 1) / (h0 - h1)
        # local minima of tau(x), periodic; parabolic sub-grid refinement
        tl, tr = np.roll(tau, 1, axis=1), np.roll(tau, -1, axis=1)
        bi, xi = np.nonzero(np.isfinite(tau) & (tau < tl) & (tau <= tr))
        t0, tm, tp = tau[bi, xi], tl[bi, xi], tr[bi, xi]
        with np.errstate(invalid='ignore', divide='ignore'):
            den = tm - 2 * t0 + tp
            ok = np.isfinite(den) & (den > 0)
            off = np.where(ok, 0.5 * (tm - tp) / den, 0.0)
            tv = np.where(ok, t0 - 0.25 * (tm - tp) * off, t0)
        out.append((bi + s, (xi + off) * dx, tv))
    return tuple(np.concatenate(v) for v in zip(*out))


# ----------------------------------------------------------------- statistics
def ring(lig, N):
    """Counts per ligament and periodic next/previous index of each site."""
    cnt = np.bincount(lig, minlength=N)
    start = np.cumsum(cnt) - cnt
    i = np.arange(lig.size)
    first, last = start[lig], start[lig] + cnt[lig] - 1
    return cnt, np.where(i == last, first, i + 1), np.where(i == first, last, i - 1)


def gaps_of(lig, xs, N):
    """Length from each site to the next one (periodic). One site alone: L."""
    cnt, nxt, _ = ring(lig, N)
    g = (xs[nxt] - xs) % L
    g[cnt[lig] == 1] = L
    return cnt, nxt, g


def lags_of(lig, taus, N):
    """tau_i minus the later of the two neighbours; NaN if fewer than 3 sites."""
    cnt, nxt, prv = ring(lig, N)
    lag = taus - np.maximum(taus[nxt], taus[prv])
    lag[cnt[lig] < 3] = np.nan
    return lag


def st(v, pcts=(5, 50, 95)):
    q = np.percentile(v, pcts)
    return (f"mean {v.mean():.3f}  std {v.std():.3f}  "
            + "  ".join(f"p{p}={a:.3f}" for p, a in zip(pcts, q)))


def dw(g):
    return (1.5 * g / W) ** (1 / 3)


def dw_peaks(g):
    """Histogram of D/W with bins of 0.05. Returns the height of the highest bin
    in 2.1 < D/W < 2.6 relative to the highest bin in 1.7 < D/W < 2.1, the
    centres of these two bins, and the histogram as fractions for 1.5..3.6."""
    hist, e = np.histogram(dw(g), bins=np.arange(0.8, 4.7001, 0.05))
    c = 0.5 * (e[1:] + e[:-1])
    w1, w2 = (c > 1.7) & (c < 2.1), (c > 2.1) & (c < 2.6)
    show = (c > 1.5) & (c < 3.6)
    return (hist[w2].max() / hist[w1].max(), c[w1][hist[w1].argmax()],
            c[w2][hist[w2].argmax()], hist[show] / g.size)


def fwhm(v, lo, hi, bw):
    hist, e = np.histogram(v, bins=np.arange(lo, hi + bw / 2, bw))
    above = np.nonzero(hist >= 0.5 * hist.max())[0]
    return e[above[-1] + 1] - e[above[0]], 0.5 * (e[hist.argmax()] + e[hist.argmax() + 1])


def items_1_to_4(lig, xs, taus, N, hist=False):
    res = {}
    first = np.full(N, np.inf)
    np.minimum.at(first, lig, taus)
    first = first[np.isfinite(first)]
    P(f"[1] first pinch time per ligament: {st(first)}   (ligaments with a pinch: {first.size}/{N})")
    cnt, nxt, g = gaps_of(lig, xs, N)
    P(f"[2] pinch sites per ligament: mean {cnt.mean():.3f}  std {cnt.std():.3f}  "
      f"min {cnt.min()}  max {cnt.max()}")
    P("    distribution of the count: " + "  ".join(
        f"{c}:{f:.3f}" for c, f in enumerate(np.bincount(cnt) / N) if f > 0.0005))
    gl = g[cnt[lig] >= 2] / LAM
    P(f"[3] gap between adjacent sites, in units of 9.02: {st(gl)}   (n = {gl.size})")
    for n in (1, 2, 3):
        P(f"    fraction of gaps within 0.25 of {n}: {np.mean(np.abs(gl-n) < 0.25):.4f}; "
          f"within 0.5 of {n}: {np.mean(np.abs(gl-n) < 0.5):.4f}")
    P(f"    fraction of gaps < 0.5: {np.mean(gl < 0.5):.4f};  > 3.5: {np.mean(gl > 3.5):.5f}")
    if hist:
        h, e = np.histogram(gl, bins=np.arange(0, 4.0001, 0.1))
        P("    histogram of gap/9.02 (bin lower edge: fraction of all gaps), bin width 0.1")
        for a in range(0, 40, 8):
            P("      " + "  ".join(f"{e[b]:.1f}:{h[b]/gl.size:.4f}" for b in range(a, a + 8)))
    main = gl[(gl > 0.5) & (gl < 1.5)]
    d = dw(main * LAM)
    w, mode = fwhm(d, 1.4, 2.3, 0.01)
    P(f"    gaps in 0.5..1.5 wavelengths ({main.size/gl.size:.4f} of all): gap/9.02 mean {main.mean():.4f} "
      f"std {main.std():.4f};  D/W mean {d.mean():.4f} std {d.std():.4f}  "
      f"p16={np.percentile(d,16):.3f} p84={np.percentile(d,84):.3f}  mode {mode:.3f}  FWHM {w:.3f}")
    lag = lags_of(lig, taus, N)
    lag = lag[np.isfinite(lag)]
    P(f"[4] lag = tau_i - later neighbour: {st(lag, (5, 25, 50, 75, 95))}   (n = {lag.size})")
    thr = (0, 0.5, 1, 2, 3, 4.5)
    fr = [np.mean(lag > t) for t in thr]
    P("    fraction of sites with lag > t:  " + "  ".join(f">{t}: {f:.4f}" for t, f in zip(thr, fr)))
    res.update(first=first.mean(), first_sd=first.std(), nsites=cnt.mean(), gap_sd=gl.std(),
               main_sd=main.std(), dw_sd=d.std(), dw_fwhm=w, lag_mean=lag.mean(),
               lag_sd=lag.std(), lag_fr=fr, near1=np.mean(np.abs(gl - 1) < 0.5))
    return res


def fragment_report(lig, xs, N, rng, label, verbose=True, nfin=None):
    """Fragment statistics for a given set of pinched sites. nfin: final mean
    number of sites, used for the reference case of necks that pinch
    independently of each other with probability p = (mean pinched)/nfin."""
    cnt, nxt, g = gaps_of(lig, xs, N)
    gl = g / LAM
    n = np.rint(gl).astype(int)
    pn = np.bincount(n, minlength=10)[:9] / n.size
    big = np.mean(n > 8)
    ratio = pn[2:9] / np.where(pn[1:8] > 0, pn[1:8], np.nan)
    multi = cnt[lig] >= 3
    if multi.sum() > 10:
        r = np.corrcoef(g[multi], g[nxt][multi])[0, 1]
        gs = g[np.lexsort((rng.random(g.size), lig))]  # shuffled inside each ligament
        r0 = np.corrcoef(gs[multi], gs[nxt][multi])[0, 1]
    else:
        r = r0 = np.nan
    P(f"{label}: pinched sites per ligament mean {cnt.mean():.3f} std {cnt.std():.3f}; "
      f"ligaments with none {np.mean(cnt == 0):.4f}; fragments {g.size}")
    P("    P(n), n=0..8: " + " ".join(f"{p:.4f}" for p in pn) + f"   P(n>8) = {big:.4f}")
    P("    P(n+1)/P(n), n=1..7: " + " ".join(f"{a:.3f}" for a in ratio))
    if nfin:
        p = cnt.mean() / nfin
        P(f"    reference, independent necks: P(1) = p = {p:.3f}, constant ratio 1-p = {1-p:.3f}")
    P(f"    correlation of successive fragment lengths {r:+.3f}; same after shuffling the "
      f"order inside each ligament {r0:+.3f}")
    if g.size > 100:
        pr, c1, c2, dh = dw_peaks(g)
        P(f"    D/W histogram (bin 0.05): highest bin in 1.7..2.1 at {c1:.3f}; highest bin in 2.1..2.6 "
          f"at {c2:.3f}, height ratio {pr:.3f};  fraction of fragments with D/W > 2.8: "
          f"{np.mean(dw(g) > 2.8):.4f};  mean length {gl.mean():.3f}")
        if verbose:
            P("    D/W histogram, bins of 0.05 from 1.5 to 3.6 (fraction of fragments): "
              + " ".join(f"{a:.3f}" for a in dh))
    if verbose:
        h, _ = np.histogram(gl, bins=np.arange(0, 8.0001, 0.25))
        P("    unrounded length/9.02, bins of 0.25 from 0 to 8: "
          + " ".join(f"{a:.3f}" for a in h / g.size))
    return cnt.mean(), pn


def by_unpinched(lig, xs, taus, N, T):
    """Fragment length at time T, grouped by the number m of necks inside the
    fragment that have not pinched yet (they pinch later in linear theory)."""
    cntf = np.bincount(lig, minlength=N)
    pos = np.arange(lig.size) - (np.cumsum(cntf) - cntf)[lig]
    s = taus <= T
    _, nxt, g = gaps_of(lig[s], xs[s], N)
    m = (pos[s][nxt] - pos[s] - 1) % cntf[lig[s]]
    P(f"T={T}: fragment length by number m of unpinched necks inside  " + ";  ".join(
        f"m={k}: fraction {np.mean(m == k):.3f}, length/9.02 mean {g[m == k].mean()/LAM:.3f} "
        f"sd {g[m == k].std()/LAM:.3f}, D/W mean {dw(g[m == k]).mean():.3f}" for k in range(4)))


# ----------------------------------------------------------------- main
def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    NC = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    rng = np.random.default_rng(SEED)
    ZA, ZB = rng.standard_normal((N, 25)), rng.standard_normal((N, 25))
    rs = np.random.default_rng(SEED + 1)  # for shuffles only
    xm, sm, s2 = dispersion_checks()

    sig = 0.08 / np.sqrt(25)
    P("=" * 78)
    P(f"CONVERGENCE CHECK on the first {NC} ligaments (base case)")
    ref = None
    for name, kw in (("production nx=2048 dt=0.02 coarse step 0.2", dict()),
                     ("brute force nx=2048 dt=0.02 every step", dict(r=1)),
                     ("refined nx=4096 dt=0.01 coarse step 0.1", dict(nx=4096, dt=0.01)),
                     ("refined nx=8192 dt=0.005 coarse step 0.05", dict(nx=8192, dt=0.005))):
        t0 = time.time()
        lg, xs, ts = find_sites(sig * ZA[:NC], sig * ZB[:NC], **kw)
        msg = f"  {name}: {lg.size} sites, {time.time()-t0:.1f} s"
        o = np.lexsort((ts, lg))  # order by pinch time inside each ligament
        lg, xs, ts = lg[o], xs[o], ts[o]
        if ref is None:
            ref = (lg, xs, ts)
        elif lg.size == ref[0].size and np.array_equal(lg, ref[0]):
            dxp = (xs - ref[1] + L / 2) % L - L / 2  # periodic distance
            msg += (f"; same number of sites in every ligament; max |dx| {np.abs(dxp).max():.2e}, "
                    f"max |dtau| {np.abs(ts-ref[2]).max():.2e}, median |dtau| "
                    f"{np.median(np.abs(ts-ref[2])):.1e}")
        else:
            msg += "; SITE LIST DIFFERS from production"
        P(msg)

    P("=" * 78)
    P(f"BASE CASE: {N} ligaments, 25 modes, rms 0.08, seed {SEED}")
    t0 = time.time()
    lig, xs, taus = find_sites(sig * ZA, sig * ZB)
    P(f"(site search took {time.time()-t0:.0f} s; rms of h(x,0) over the sample: "
      f"{sig*np.sqrt((ZA**2+ZB**2).sum(1).mean()/2):.4f})")
    base = items_1_to_4(lig, xs, taus, N, hist=True)
    ts = np.sort(taus)
    nbar = lambda T: np.searchsorted(ts, T, side='right') / N
    P("    mean number of sites with tau <= T: " + "  ".join(
        f"T={T}: {nbar(T):.3f}" for T in (15, 20, 25, 30, 35, 40)))
    lag = lags_of(lig, taus, N)

    P("-" * 78)
    P("[5] TIME-RESOLVED FRAGMENTS (n = length/9.02 rounded; a ligament with one pinched "
      "site is one fragment of length 100)")
    for T in (8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 40):
        s = taus <= T
        fragment_report(lig[s], xs[s], N, rs, f"T={T}", nfin=base['nsites'])

    for T in (11, 12, 14):
        by_unpinched(lig, xs, taus, N, T)

    P("-" * 78)
    P("[6] COMPARISON WITH THE SIMULATION")
    nfin = base['nsites']
    for target in (11.73 / 2, 13.87 / 2):
        Tm = ts[int(target * N)]
        P(f"mean number of pinched sites = {target:.3f} at linear-theory time T = {Tm:.3f}")
        s = taus <= Tm
        fragment_report(lig[s], xs[s], N, rs, f"   linear theory at T={Tm:.2f}", nfin=nfin)
    P(f"mean pinched sites at T=12: {nbar(12):.3f}, at T=14: {nbar(14):.3f}; final (t<=40): {nfin:.3f}")
    lf = lag[np.isfinite(lag)]
    for target in (13.87 / 2, 11.73 / 2):
        f = 1 - target / nfin
        th = np.quantile(lf, 1 - f)
        P(f"to end with {target:.3f} fragments, a fraction {f:.4f} of the {nfin:.3f} necks must "
          f"fail; the fraction of sites with lag > {th:.3f} equals that")
        keep = ~(lag > th)
        fragment_report(lig[keep], xs[keep], N, rs,
                        f"   all necks pinch except those with lag > {th:.2f}", verbose=False)
    for th in (0.0, 1.0, 2.0):
        keep = ~(lag > th)
        fragment_report(lig[keep], xs[keep], N, rs,
                        f"   all necks pinch except those with lag > {th:.1f}", verbose=False)
    P("Time-resolved test of the neck-failure rule: sites with tau <= T, without the "
      "necks that have lag > threshold")
    for th in (0.0, -0.5):
        keep = ~(lag > th)
        P(f"  threshold {th}: mean number of pinched sites  " + "  ".join(
            f"T={T}: {np.sum(keep & (taus <= T))/N:.3f}" for T in (10, 11, 12, 13, 14, 16, 18, 40)))
        for T in (12, 14, 18):
            s = keep & (taus <= T)
            fragment_report(lig[s], xs[s], N, rs, f"   threshold {th}, T={T}", verbose=(th == 0.0))

    P("=" * 78)
    P("[7] SENSITIVITY (same random numbers as the base case, rescaled)")
    n = np.arange(1, 26)
    cases = [
        ("(a) rms 0.02, 25 modes", np.full(25, 0.02 / 5), 48.0),
        ("(a) rms 0.005, 25 modes", np.full(25, 0.005 / 5), 56.0),
        ("(b) 15 unstable modes, same variance per mode as base (rms 0.062)", np.full(15, sig), 40.0),
        ("(b') 15 unstable modes, total rms 0.08", np.full(15, 0.08 / np.sqrt(15)), 40.0),
        ("(c) variance ~ n^-2, same variance at n=11 as base", sig * 11 / n, 40.0),
    ]
    summ = [("base: rms 0.08, 25 modes", base)]
    for name, sv, tmax in cases:
        P("-" * 78)
        nm = sv.size
        P(f"{name}; rms of h(x,0) = {np.sqrt((sv**2).sum()):.4f}; t from 0 to {tmax}")
        t0 = time.time()
        lg, x_, t_ = find_sites(sv * ZA[:, :nm], sv * ZB[:, :nm], tmax=tmax)
        P(f"(site search took {time.time()-t0:.0f} s)")
        summ.append((name, items_1_to_4(lg, x_, t_, N)))
    P("-" * 78)
    P("Summary. sd = standard deviation; gaps in units of 9.02; 'main' = gaps in 0.5..1.5.")
    P("bandwidth = 1/sqrt(2 |s''| t1)/x_max: relative spectral width of the amplified "
      "noise at the mean first-pinch time t1")
    P(f"{'case':66s} {'t1':>6s} {'sites':>6s} {'gap sd':>7s} {'main sd':>7s} {'D/W sd':>7s} "
      f"{'bandw':>6s} {'lag mean':>8s} {'lag sd':>6s} {'lag>0':>6s} {'lag>1':>6s} {'lag>2':>6s}")
    b0 = summ[0][1]
    for name, r in summ[1:3]:
        e = np.log(r['first'] / b0['first'])
        P(f"exponent p in (width ~ t1^p) between base and {name}: gap sd {np.log(r['gap_sd']/b0['gap_sd'])/e:+.3f}, "
          f"main sd {np.log(r['main_sd']/b0['main_sd'])/e:+.3f}, D/W sd {np.log(r['dw_sd']/b0['dw_sd'])/e:+.3f} "
          f"(bandwidth estimate: -0.5)")
    for name, r in summ:
        bw = 1 / np.sqrt(2 * abs(s2) * r['first']) / xm
        P(f"{name:66s} {r['first']:6.2f} {r['nsites']:6.2f} {r['gap_sd']:7.4f} {r['main_sd']:7.4f} "
          f"{r['dw_sd']:7.4f} {bw:6.3f} {r['lag_mean']:8.3f} {r['lag_sd']:6.3f} {r['lag_fr'][0]:6.3f} "
          f"{r['lag_fr'][2]:6.3f} {r['lag_fr'][3]:6.3f}")


if __name__ == "__main__":
    main()
