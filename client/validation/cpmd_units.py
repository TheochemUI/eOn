# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy>=1.12"]
# ///
"""Exact checks of the unit chain between eOn and CPMD.

Run with:  uv run --with sympy client/validation/cpmd_units.py

eOn works in eV and Angstrom. CPMD reads Angstrom (ANGSTROM keyword) and
converts with its own Bohr, cnst.mod fbohr = 1/0.529177210859, then answers
in Hartree and Hartree/Bohr. Both routes to CPMD convert back with CODATA
2018 (27.211386245988 eV, 0.529177210903 A): cpmd_extpot.py for the file
route, rgpot units.hpp and cpmdc for the in-process route.

1. The force conversion. The exact force is dE/dx_A = (dE/dx_b) / a0_cpmd,
   because CPMD's x_b is x_A / a0_cpmd. Using a0_2018 instead leaves a
   relative error a0_cpmd / a0_2018 - 1, printed here (-8.3e-11), far
   below the SCF tolerance of a force.
2. The sign. CPMD's fion (GEOMETRY columns 4-6, the "GRADIENTS (-FORCES)"
   table) is the force -dE/dx. cpmdc hands rgpot the gradient -fion, and
   rgpot multiplies the gradient by NEG_GRAD_TO_FORCE = -Ha/a0, so both
   routes return +fion * Ha/a0. Checked on E = k |x|^2.
3. The stress. cpmdc returns paiu/omega in Ha/Bohr^3. CPMD prints
   (paiu/omega) * au_kb with au_kb = 294210.1080 kbar per Ha/Bohr^3
   (cnst.mod); CODATA 2018 gives 294210.157, a relative difference of
   -1.7e-7. rgpot's HARTREE_PER_BOHR3_TO_EV_PER_ANGSTROM3 is Ha/a0^3.
4. The strain derivative. sigma_xx = (1/V) dE/d eps_xx at fixed
   fractional coordinates; the central difference over +-eps carries the
   error eps^2/6 E'''(0)/V: 4.2e-6 E'''/V at +-0.005.
"""

import sympy as sp

R = sp.Rational
HA_EV = R("27.211386245988")
A0_2018 = R("0.529177210903")
A0_CPMD = R("0.529177210859")
AU_KB_CPMD = R("294210.1080")
# 1 Ha/Bohr^3 in Pa: E_h / a0^3 with E_h = 4.3597447222071e-18 J.
EH_J = R("4.3597447222071e-18")


def check_force_factor():
    x_a, k = sp.symbols("x_A k", positive=True)
    x_b = x_a / A0_CPMD
    energy_ha = k * x_b**2
    exact = sp.diff(-energy_ha * HA_EV, x_a)
    fion = sp.diff(-energy_ha, x_a) * A0_CPMD  # -dE/dx_b, Ha/Bohr_cpmd
    used = fion * HA_EV / A0_2018
    rel = sp.simplify(used / exact - 1)
    assert rel == A0_CPMD / A0_2018 - 1
    print(f"force factor used        {sp.N(HA_EV / A0_2018, 15)} eV/A per Ha/Bohr")
    print(f"relative error (a0 pair) {sp.N(rel, 3)}")
    assert abs(rel) < R(1, 10**10)


def check_sign():
    x, k = sp.symbols("x k", positive=True)
    energy = k * x**2
    fion = -sp.diff(energy, x)          # CPMD force
    grad_cpmdc = -fion                  # cpmd_embed_c_api: -fion
    neg_grad_to_force = -HA_EV / A0_2018
    rgpot_force = grad_cpmdc * neg_grad_to_force
    extpot_force = fion * HA_EV / A0_2018
    assert sp.simplify(rgpot_force - extpot_force) == 0
    assert sp.simplify(rgpot_force + sp.diff(energy, x) * HA_EV / A0_2018) == 0
    print("sign: both routes return +fion * Ha/a0 = -dE/dx")


def check_stress():
    au_kb_2018 = EH_J / (A0_2018 * R(1, 10**10)) ** 3 / R(10**8)
    rel = AU_KB_CPMD / au_kb_2018 - 1
    print(f"au_kb CPMD {AU_KB_CPMD}, CODATA 2018 {sp.N(au_kb_2018, 12)}, rel {sp.N(rel, 3)}")
    assert abs(rel) < R(1, 10**6)
    ev_a3 = HA_EV / A0_2018**3
    print(f"Ha/Bohr^3 to eV/A^3      {sp.N(ev_a3, 15)}")


def check_strain_difference():
    eps, e0, e1, e2, e3, V = sp.symbols("epsilon E0 E1 E2 E3 V", real=True)
    s = sp.Symbol("s", real=True)
    E = e0 + e1 * s + e2 * s**2 / 2 + e3 * s**3 / 6
    central = (E.subs(s, eps) - E.subs(s, -eps)) / (2 * eps)
    err = sp.expand(central - sp.diff(E, s).subs(s, 0))
    assert sp.simplify(err - e3 * eps**2 / 6) == 0
    print("central strain difference error: eps^2/6 E'''")


if __name__ == "__main__":
    check_force_factor()
    check_sign()
    check_stress()
    check_strain_difference()
