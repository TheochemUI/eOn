# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy>=1.12"]
# ///
"""Symbolic checks for the formulas the improved dimer and the
finite-difference Hessian-vector products use.

Run with:  uv run --with sympy client/validation/dimer_rotation.py

1. The forward-image gradient at the optimal rotation angle, interpolated
   from g0, g1 and the trial gradient g1' (ImprovedDimer::compute, "saves
   one force call"), is exact on a quadratic surface.
2. The two-term Fourier model of the curvature, with a1 and b1 taken from
   one trial rotation, is exact on a quadratic surface, and phi_min is its
   stationary point.
3. The forward-difference curvature carries an O(delta) error led by the
   third derivative; the central and fourth-order force derivatives in
   FiniteDifference.h have O(h^2) and O(h^4) leading errors.
"""

import sympy as sp


def check_interpolation():
    phi, phip, d = sp.symbols("phi phi_p delta", real=True)
    n = 3
    H = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"h{min(i, j)}{max(i, j)}"))
    g0 = sp.Matrix(sp.symbols("g0_0:3"))
    # Orthonormal tau, theta in the rotation plane.
    tau = sp.Matrix([1, 0, 0])
    theta = sp.Matrix([0, 1, 0])

    def g1(angle):
        # Gradient at x0 + delta tau(angle) on E = g0.x + x.H.x / 2.
        return g0 + d * H * (tau * sp.cos(angle) + theta * sp.sin(angle))

    g1_0, g1_p = g1(0), g1(phip)
    # ImprovedDimer.cpp interpolation, term by term.
    interp = (
        g1_0 * sp.sin(phip - phi) / sp.sin(phip)
        + g1_p * sp.sin(phi) / sp.sin(phip)
        + g0 * (1 - sp.cos(phi) - sp.sin(phi) * sp.tan(phip / 2))
    )
    # Every trig term rewritten in t = tan(phi_p / 2) and u = tan(phi / 2)
    # makes the residual rational, so cancel() decides it exactly.
    t, u = sp.symbols("t u", real=True)
    weierstrass = {
        sp.sin(phi): 2 * u / (1 + u**2),
        sp.cos(phi): (1 - u**2) / (1 + u**2),
        sp.sin(phip): 2 * t / (1 + t**2),
        sp.cos(phip): (1 - t**2) / (1 + t**2),
        sp.tan(phip / 2): t,
    }
    diff = (interp - g1(phi)).applyfunc(
        lambda e: sp.cancel(sp.expand_trig(e).xreplace(weierstrass))
    )
    assert all(e == 0 for e in diff), diff
    print("g1(phi_min) interpolation: exact on a quadratic surface")


def check_fourier_curvature():
    phi, phip = sp.symbols("phi phi_p", real=True)
    a, b, c = sp.symbols("a b c", real=True)  # 2x2 Hessian in the plane
    H = sp.Matrix([[a, b], [b, c]])

    def C(angle):
        t = sp.Matrix([sp.cos(angle), sp.sin(angle)])
        return (t.T * H * t)[0]

    C0, Cp = C(0), C(phip)
    dC0 = sp.diff(C(phi), phi).subs(phi, 0)
    b1 = dC0 / 2
    a1 = (C0 - Cp + b1 * sp.sin(2 * phip)) / (1 - sp.cos(2 * phip))
    a0 = 2 * (C0 - a1)
    model = a0 / 2 + a1 * sp.cos(2 * phi) + b1 * sp.sin(2 * phi)
    resid = sp.simplify(sp.expand_trig(model - C(phi)))
    assert resid == 0, resid
    # phi_min = atan(b1/a1)/2 is stationary for the model.
    phimin = sp.atan(b1 / a1) / 2
    stat = sp.simplify(sp.diff(model, phi).subs(phi, phimin))
    assert stat == 0, stat
    # The code's dC/dphi at 0 is 2 (g1 - g0).theta / delta = 2 theta.H.tau.
    assert sp.simplify(dC0 - 2 * b) == 0
    print("Fourier curvature model: exact, phi_min stationary")


def check_fd_orders():
    s, h, d = sp.symbols("s h delta", positive=True)
    # Energy as its Taylor polynomial about x: e[k] is the k-th derivative.
    e = sp.symbols("e0:9")
    E = sum(e[k] * s**k / sp.factorial(k) for k in range(9))
    f = -sp.diff(E, s)  # force along the step

    def at(step):
        return sp.expand(f.subs(s, step))

    # Forward-difference curvature (g(x+d) - g(x)) / d = -(f(d) - f(0)) / d.
    fwd = sp.expand(-(at(d) - at(0)) / d)
    lead = sp.series(fwd - e[2], d, 0, 2).removeO()
    assert sp.expand(lead - d / 2 * e[3]) == 0, lead
    exact = sp.diff(f, s).subs(s, 0)  # f'(0) = -e2
    central = (at(h) - at(-h)) / (2 * h)
    fourth = (-at(2 * h) + 8 * at(h) - 8 * at(-h) + at(-2 * h)) / (12 * h)
    ec = sp.series(sp.expand(central - exact), h, 0, 4).removeO()
    e4 = sp.series(sp.expand(fourth - exact), h, 0, 6).removeO()
    # The force's third and fifth derivatives are -e4 and -e6.
    assert sp.expand(ec - h**2 / 6 * (-e[4])) == 0, ec
    assert sp.expand(e4 + h**4 / 30 * (-e[6])) == 0, e4
    print("forward curvature error delta/2 E3; central h^2/6 f3; fourth -h^4/30 f5")


if __name__ == "__main__":
    check_interpolation()
    check_fourier_curvature()
    check_fd_orders()
