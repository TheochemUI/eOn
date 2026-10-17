# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy>=1.12"]
# ///
"""Fixed points of the two NEB parallel spring forces, and the effect of
removing the net translation from the projected force only.

Run with:  uv run --with sympy client/validation/neb_spring_fixed_point.py

The spring forces as the code writes them, with d+ = R[i+1] - R[i],
d- = R[i] - R[i-1] and the unit tangent tau:

- Hooke (UniformSpring, the default since 6a89efd32):
  F_par = k (|d+| - |d-|) tau, the form of Jonsson, Mills and Jacobsen
  (1998) and of Henkelman, Uberuaga and Jonsson (2000).
- squared (GeometricSpring, the SVN form at svn-Jul_01_2024, where
  distPrev and distNext are squaredNorm()):
  F_par = k (|d+|^2 - |d-|^2) tau.

1. Both vanish exactly when |d+| = |d-|; the squared form is the Hooke
   form times |d+| + |d-|, so its stiffness at the fixed point is
   2 k |d| and falls with the image spacing.
2. On a curved path (a circular arc) equal chords are equal arc lengths,
   so both laws place the converged images at equal arc length.
3. With every atom free the client removes the net translation from the
   projected force. If the tangent itself carries a translation, the
   projected force Pi (F_perp + s tau) vanishes with s != 0: the band
   "converges" with unequal spacing, and the spacing depends on the spring
   law. With the tangent and the distances made translation-free first,
   the projected force vanishes only for F_perp = 0 and s = 0.
"""

import sympy as sp


def check_fixed_points():
    k, a, b = sp.symbols("k a b", positive=True)
    hooke = k * (a - b)
    squared = k * (a**2 - b**2)
    assert sp.solve(sp.Eq(hooke, 0), a) == [b]
    assert sp.solve(sp.Eq(squared, 0), a) == [b]
    assert sp.simplify(squared - hooke * (a + b)) == 0
    stiff_h = sp.diff(hooke, a).subs(a, b)
    stiff_s = sp.diff(squared, a).subs(a, b)
    assert stiff_h == k and sp.simplify(stiff_s - 2 * k * b) == 0
    print("1. Hooke and squared springs both vanish only at |d+| = |d-|;")
    print(f"   stiffness at the fixed point: Hooke {stiff_h}, squared {stiff_s}")


def check_curved_path():
    r, t = sp.symbols("r theta", positive=True)
    chord = 2 * r * sp.sin(t / 2)
    slope = sp.diff(chord, t)
    # The chord grows strictly with the subtended angle on 0 < theta < pi,
    # so equal chords subtend equal angles and so equal arcs r*theta.
    assert sp.simplify(slope - r * sp.cos(t / 2)) == 0
    assert all(
        slope.subs({r: 1, t: v}) > 0 for v in (sp.pi / 100, sp.pi / 2, 99 * sp.pi / 100)
    )
    print("2. on a circular arc the chord 2 r sin(theta/2) rises strictly with theta")
    print("   on (0, pi): equal chords are equal arcs")


def check_translation():
    f, s = sp.symbols("f s", real=True)
    # Two atoms in one dimension. Only atom 0 moves between the images, so
    # tau carries a translation; the true force has no net translation.
    tau = sp.Matrix([1, 0])
    F = sp.Matrix([f, -f])
    Pi = sp.eye(2) - sp.ones(2, 2) / 2
    fperp = F - (F.dot(tau)) * tau
    projected = Pi * (fperp + s * tau)
    sol = sp.solve(list(projected), s, dict=True)
    assert sol == [{s: -f}]
    print("3. translation removed from the projected force only:")
    print(f"   Pi (F_perp + s tau) = 0 at s = {sol[0][s]}, a nonzero spring force")
    taup = Pi * tau
    taup = taup / sp.sqrt(taup.dot(taup))
    fperp2 = F - (F.dot(taup)) * taup
    assert sp.simplify(fperp2) == sp.zeros(2, 1)
    projected2 = Pi * (fperp2 + s * taup)
    sol2 = sp.solve(list(projected2), s, dict=True)
    assert sol2 == [{s: 0}]
    print("   with a translation-free tangent: F_perp = 0 and s = 0 only")

    # The same on 13 atoms in three dimensions with rational data.
    n = 13
    rng = sp.Rational
    tau13 = sp.Matrix([rng((7 * i) % 11 - 5, 3) for i in range(3 * n)])
    tau13[0] += 4  # one atom moves much more: a net translation
    F13 = sp.Matrix([rng((5 * i) % 13 - 6, 2) for i in range(3 * n)])
    # The true force of an isolated cluster has no net translation.
    for c in range(3):
        mean = sum(F13[3 * j + c] for j in range(n)) / n
        for j in range(n):
            F13[3 * j + c] -= mean

    def remove_translation(v):
        w = v.copy()
        for c in range(3):
            mean = sum(w[3 * j + c] for j in range(n)) / n
            for j in range(n):
                w[3 * j + c] -= mean
        return w

    def residual(t):
        t = t / sp.sqrt(t.dot(t))
        fp = F13 - F13.dot(t) * t
        # Least-squares s for the remaining projected force.
        a = remove_translation(t)
        b = remove_translation(fp)
        s_best = -b.dot(a) / a.dot(a)
        return s_best, remove_translation(fp + s_best * t)

    s_raw, _ = residual(tau13)
    s_free, _ = residual(remove_translation(tau13))
    assert s_raw != 0 and s_free == 0
    print(f"   13 atoms: best spring force {sp.N(s_raw, 6)} with the raw tangent,")
    print(f"   {s_free} with the translation-free tangent")


if __name__ == "__main__":
    check_fixed_points()
    check_curved_path()
    check_translation()
