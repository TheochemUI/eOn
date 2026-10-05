# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy>=1.12", "mpmath>=1.3", "numpy>=1.26"]
# ///
"""Reference values and symbolic checks for the ring-polymer instanton,
the path-integral ring and its thermostat, independent of eOn.

Run with:  uv run client/validation/instanton_references.py

Units are eOn's: eV, amu, Angstrom, kelvin; hbar = 0.06465415129579072
eV^0.5 amu^0.5 Angstrom, one time unit sqrt(amu Angstrom^2 / eV).

Symbolic (sympy):
  1. The ring potential U_N, its gradient and its Hessian blocks
     (H_j + 2c on the diagonal, -c beside it and in the closing corner).
  2. The half-ring folding: on rings with q_j = q_{N-j}, U_N / 2 written
     over beads 0..N/2 has the gradient and the chain blocks the half
     ring uses (H/2 + c at the two turning points, -c off the diagonal).
  3. The free-ring normal-mode frequencies 2 omega_N sin(pi k / N).
  4. The PILE step keeps the Maxwell-Boltzmann variance N kB T m and has
     the critical damping time 1 / (2 omega_k).
  5. The splitting zero-mode identity det' J = det J (v^T J^-1 v).

Numerical references (mpmath at 40 digits, numpy for the DVR):
  A. Symmetric Eckart V0 sech^2(x / a), V0 = 0.425 eV, a = 0.734 amu^0.5
     Angstrom, unit mass: T_c = hbar sqrt(2 V0) / (2 pi kB a), the exact
     flux (1 / 2 pi hbar) int P(E) exp(-beta E) dE and the N -> infinity
     instanton, the steepest-descent value of the same integral with
     P = exp(-theta(E)), theta = 2 sqrt(2) pi a (sqrt(V0) - sqrt(E)) / hbar.
  B. The harmonic reactant: ln Z_N = -sum_k ln(beta_N hbar omega_k) with
     omega_k^2 = omega^2 + 4 sin^2(pi k / N) / (beta_N hbar)^2, against
     ln Z = -ln(2 sinh(beta hbar omega / 2)); the error tends to
     u^3 coth(u / 2) / (48 N^2), u = beta hbar omega.
  C. The quartic double well V0 (x^2 - 1)^2, unit mass: the exact
     splitting from a sinc DVR (Colbert and Miller) and the continuum
     instanton 2 hbar omega sqrt(4 omega / (pi hbar)) exp(-S0 / hbar),
     S0 = 2 omega / 3, omega^2 = 8 V0.
  D. Round-off: ln(2 sinh(x / 2)) as eOn writes it and the parabolic
     factor near T_c, against mpmath.
"""

import math

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 40
HBAR = mp.mpf("0.06465415129579072")
KB = mp.mpf("8.617333262145177e-5")


# ---------------------------------------------------------------- symbolic
def check_ring_derivatives(n=6):
    q = sp.symbols(f"q0:{n}", real=True)
    c = sp.symbols("c", positive=True)
    V = sp.Function("V")
    u = sum(V(qj) for qj in q) + c / 2 * sum(
        (q[(j + 1) % n] - q[j]) ** 2 for j in range(n)
    )
    for j in range(n):
        g = sp.diff(u, q[j])
        want = sp.diff(V(q[j]), q[j]) + c * (2 * q[j] - q[j - 1] - q[(j + 1) % n])
        assert sp.simplify(g - want) == 0
        for k in range(n):
            h = sp.diff(u, q[j], q[k])
            if k == j:
                assert sp.simplify(h - sp.diff(V(q[j]), q[j], 2) - 2 * c) == 0
            elif (k - j) % n in (1, n - 1):
                assert sp.simplify(h + c) == 0
            else:
                assert h == 0
    print(f"1. U_N gradient and Hessian blocks hold for N = {n}")


def check_half_ring(n=8):
    half = n // 2
    p = sp.symbols(f"p0:{half + 1}", real=True)
    c = sp.symbols("c", positive=True)
    V = sp.Function("V")
    full = [p[j] if j <= half else p[n - j] for j in range(n)]
    u = sum(V(x) for x in full) + c / 2 * sum(
        (full[(j + 1) % n] - full[j]) ** 2 for j in range(n)
    )
    w = u / 2
    for j in range(half + 1):
        g = sp.diff(w, p[j])
        turning = j in (0, half)
        vp = sp.diff(V(p[j]), p[j])
        if turning:
            nb = p[1] if j == 0 else p[half - 1]
            want = vp / 2 + c * (p[j] - nb)
        else:
            want = vp + c * (2 * p[j] - p[j - 1] - p[j + 1])
        assert sp.simplify(g - want) == 0
        for k in range(half + 1):
            h = sp.diff(w, p[j], p[k])
            vpp = sp.diff(V(p[j]), p[j], 2)
            if k == j:
                want = vpp / 2 + c if turning else vpp + 2 * c
                assert sp.simplify(h - want) == 0
            elif abs(k - j) == 1:
                assert sp.simplify(h + c) == 0
            else:
                assert h == 0
    print(f"2. half ring of N = {n}: turning blocks H/2 + c, interior H + 2c, "
          "off-diagonal -c")


def check_free_ring(n=6):
    c = sp.symbols("c", positive=True)
    t = sp.zeros(n, n)
    for j in range(n):
        t[j, j] = 2 * c
        t[j, (j + 1) % n] -= c
        t[j, (j - 1) % n] -= c
    for k in range(n):
        v = sp.Matrix([sp.cos(2 * sp.pi * k * j / n) for j in range(n)])
        lam = 4 * c * sp.sin(sp.pi * k / n) ** 2
        assert sp.simplify((t * v - lam * v).applyfunc(sp.expand_trig)) == \
            sp.zeros(n, 1)
    print(f"3. free ring of N = {n}: omega_k^2 = 4 omega_N^2 sin^2(pi k / N), "
          "omega_N^2 = c = 1 / (beta_N hbar)^2")


def check_pile():
    h, tau, kt, m = sp.symbols("h tau kT_N m", positive=True)
    d = sp.exp(-h / tau)
    noise2 = kt * (1 - d**2)  # on p / sqrt(m), as eOn draws it
    var = sp.symbols("s2", positive=True)
    # One step p' = d p + noise xi keeps var = kT_N fixed.
    fixed = sp.solve(sp.Eq(var, d**2 * var + noise2), var)[0]
    assert sp.simplify(fixed - kt) == 0
    # Critical damping of a harmonic mode: gamma = 2 omega, tau = 1 / gamma.
    om, g, s = sp.symbols("omega gamma s", positive=True)
    roots = sp.solve(s**2 + g * s + om**2, s)
    crit = sp.solve(sp.Eq(roots[0], roots[1]), g)
    assert crit == [2 * om]
    print("4. PILE: OU step keeps p^2 / m = N kB T; critical gamma = 2 omega_k, "
          "tau_k = 1 / (2 lambda omega_k)")


def check_zero_mode():
    a, b, lam0 = sp.symbols("a b lambda0", real=True)
    j = sp.diag(lam0, a, b)
    v = sp.Matrix([1, 0, 0])
    assert sp.simplify(j.det() * (v.T * j.inv() * v)[0] - a * b) == 0
    print("5. det' J = det J (v^T J^-1 v) for v the eigenvector left out")


# ---------------------------------------------------------------- numerical
def eckart_references():
    v0, a = mp.mpf("0.425"), mp.mpf("0.734")
    omega_b = mp.sqrt(2 * v0) / a
    tc = HBAR * omega_b / (2 * mp.pi * KB)
    print(f"\nA. Eckart: omega_b = {mp.nstr(omega_b, 17)}, "
          f"T_c = {mp.nstr(tc, 17)} K")
    d2 = 8 * v0 * a**2 / HBAR**2 - 1

    def p_exact(e):
        alpha = a * mp.sqrt(2 * e) / HBAR
        ca = mp.cosh(2 * mp.pi * alpha)
        return (ca - 1) / (ca + mp.cosh(mp.pi * mp.sqrt(d2)))

    def theta(e):
        return 2 * mp.sqrt(2) * mp.pi * a * (mp.sqrt(v0) - mp.sqrt(e)) / HBAR

    # theta against the quadrature it closes, at one energy.
    e = v0 / 3
    x0 = a * mp.acosh(mp.sqrt(v0 / e))
    quad = 2 / HBAR * mp.quad(
        lambda x: mp.sqrt(2 * (v0 / mp.cosh(x / a) ** 2 - e)), [-x0, 0, x0])
    assert abs(quad - theta(e)) < mp.mpf("1e-25")
    rows = []
    for frac in ("0.5", "0.35", "0.7", "0.9"):
        t = mp.mpf(frac) * tc
        beta = 1 / (KB * t)
        exact = mp.log(mp.quad(lambda en: p_exact(en) * mp.exp(-beta * en),
                               [0, v0 / 4, v0 / 2, v0, v0 + 40 / beta])
                       / (2 * mp.pi * HBAR))
        es = 2 * mp.pi**2 * a**2 / (HBAR * beta) ** 2
        th2 = mp.sqrt(2) * mp.pi * a / (2 * HBAR * es ** mp.mpf(1.5))
        inst = (-mp.log(2 * mp.pi * HBAR) - beta * es - theta(es)
                + mp.log(mp.sqrt(2 * mp.pi / th2)))
        rows.append((frac, exact, inst))
        print(f"   T = {frac} T_c: ln(k Z_r) exact {mp.nstr(exact, 15)}, "
              f"instanton (N -> inf) {mp.nstr(inst, 15)}, "
              f"ratio {mp.nstr(mp.exp(inst - exact), 8)}, "
              f"E* = {mp.nstr(es, 10)} eV")
    return rows


def harmonic_references(omega=1.0, temperature=50.0):
    w = mp.mpf(omega)
    beta = 1 / (KB * temperature)
    exact = -mp.log(2 * mp.sinh(beta * HBAR * w / 2))
    print(f"\nB. Harmonic reactant, omega = {omega}, T = {temperature} K, "
          f"beta hbar omega = {mp.nstr(beta * HBAR * w, 12)}: "
          f"ln Z = {mp.nstr(exact, 17)}")
    prev = None
    for n in (8, 16, 32, 64, 128, 256):
        bn = beta / n
        s = mp.mpf(0)
        for k in range(n):
            wk2 = w**2 + 4 * mp.sin(mp.pi * k / n) ** 2 / (bn * HBAR) ** 2
            s -= mp.log(bn * HBAR * mp.sqrt(wk2))
        err = s - exact
        order = "" if prev is None else f", order {mp.nstr(mp.log(prev / err, 2), 6)}"
        print(f"   N = {n:4d}: ln Z_N = {mp.nstr(s, 17)}, "
              f"error {mp.nstr(err, 6)}, N^2 error {mp.nstr(err * n * n, 8)}"
              f"{order}")
        prev = err
    u = beta * HBAR * w
    print(f"   leading N^2 error u^3 coth(u / 2) / 48, u = beta hbar omega: "
          f"{mp.nstr(u**3 / mp.tanh(u / 2) / 48, 8)}")


def dvr_splitting(v0, n, half=2.4):
    x = np.linspace(-half, half, n)
    dx = x[1] - x[0]
    hb = float(HBAR)
    i = np.arange(n)
    d = i[:, None] - i[None, :]
    with np.errstate(divide="ignore"):
        t = np.where(d == 0, np.pi**2 / 3, 2.0 * (-1.0) ** d / d.astype(float) ** 2)
    h = hb * hb / (2 * dx * dx) * t + np.diag(v0 * (x * x - 1) ** 2)
    e = np.linalg.eigvalsh(h)
    return e[1] - e[0]


def double_well_references():
    print("\nC. Quartic double well V0 (x^2 - 1)^2, unit mass")
    out = {}
    for v0 in (0.12, 0.2, 0.3):
        exact = [dvr_splitting(v0, n) for n in (301, 401, 601)]
        w = mp.sqrt(8 * mp.mpf(v0))
        s0 = 2 * w / 3
        inst = 2 * HBAR * w * mp.sqrt(4 * w / (mp.pi * HBAR)) * mp.exp(-s0 / HBAR)
        out[v0] = (exact[-1], inst)
        print(f"   V0 = {v0}: DVR {exact[0]:.10e} (301) {exact[1]:.10e} (401) "
              f"{exact[2]:.10e} (601); instanton (P -> inf) "
              f"{mp.nstr(inst, 12)}, ratio {mp.nstr(inst / exact[-1], 8)}, "
              f"S0 / hbar = {mp.nstr(s0 / HBAR, 8)}")
    return out


def roundoff():
    print("\nD. Round-off")
    worst = 0.0
    for x in (1e-12, 1e-9, 1e-6, 1e-3, 0.1, 1.0, 30.0, 700.0):
        mine = 0.5 * x + math.log1p(-math.exp(-x))
        ref = float(mp.log(2 * mp.sinh(mp.mpf(x) / 2)))
        worst = max(worst, abs(mine - ref))
        print(f"   ln 2 sinh(x/2), x = {x:g}: eOn {mine:.17g}, "
              f"exact {ref:.17g}, |diff| {abs(mine - ref):.2e}")
    for rel in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10):
        tc, t = 150.0, 150.0 * (1 + rel)
        phase = math.pi * tc / t
        mine = phase / math.sin(phase)
        p = mp.pi * mp.mpf(tc) / mp.mpf(t)
        ref = p / mp.sin(p)
        print(f"   parabolic factor, T / T_c - 1 = {rel:g}: eOn {mine:.17g}, "
              f"relative error {float(abs(mine - ref) / ref):.2e}")
    return worst


if __name__ == "__main__":
    check_ring_derivatives()
    check_half_ring()
    check_free_ring()
    check_pile()
    check_zero_mode()
    eckart_references()
    harmonic_references()
    double_well_references()
    roundoff()
