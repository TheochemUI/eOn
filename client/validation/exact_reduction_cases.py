# /// script
# requires-python = ">=3.10"
# ///
"""Reference cases for eonc::exact::dot and eonc::exact::sum.

Run with:  python3 client/validation/exact_reduction_cases.py > \\
    client/unit_tests/data/systems/exact_reduction/cases.txt

Each expected value is the exact rational sum (fractions.Fraction) rounded
once to the nearest double, ties to even: CPython's int / int true
division is correctly rounded, so float(Fraction) is the reference
rounding. A sum past the overflow threshold is written as inf.

File format, one case per block, every double as a C99 hex float:

    dot <n>        or    sum <n>
    <x_0> ... <x_{n-1}>
    <y_0> ... <y_{n-1}>   (dot only)
    <expected>

The cases cover the ways a floating-point reduction loses information:
catastrophic cancellation, operands 2^1000 apart, subnormal products,
ties at the last bit (with and without a sticky bit below), sums at the
overflow threshold, and random vectors of the shapes the optimizers use.
"""

from __future__ import annotations

import random
import sys
from fractions import Fraction

MAX = float.fromhex("0x1.fffffffffffffp+1023")
MIN_SUB = float.fromhex("0x1p-1074")


def rn(q: Fraction) -> float:
    try:
        return float(q)
    except OverflowError:
        return float("inf") if q > 0 else float("-inf")


def exact_dot(x: list[float], y: list[float]) -> float:
    return rn(sum((Fraction(a) * Fraction(b) for a, b in zip(x, y)), Fraction(0)))


def exact_sum(x: list[float]) -> float:
    return rn(sum((Fraction(a) for a in x), Fraction(0)))


def emit_dot(x: list[float], y: list[float]) -> None:
    print(f"dot {len(x)}")
    print(" ".join(v.hex() for v in x))
    print(" ".join(v.hex() for v in y))
    print(exact_dot(x, y).hex())


def emit_sum(x: list[float]) -> None:
    print(f"sum {len(x)}")
    print(" ".join(v.hex() for v in x))
    print(exact_sum(x).hex())


def main() -> None:
    rng = random.Random(20261018)
    one = 1.0
    ulp = 2.0**-52

    # Cancellation: the plain left-to-right sum is 0, the exact sum is 1.
    emit_sum([1e100, 1.0, -1e100])
    emit_dot([1e100, 1.0, -1e100], [1.0, 1.0, 1.0])
    # Ties at the last bit: to even, then broken by a sticky bit far below.
    emit_sum([one, ulp / 2])
    emit_sum([one + ulp, ulp / 2])
    emit_sum([one, ulp / 2, 2.0**-200])
    emit_sum([one, ulp / 2, -(2.0**-200)])
    # Subnormal results and products that underflow on their own.
    emit_sum([MIN_SUB, MIN_SUB, MIN_SUB])
    emit_dot([2.0**-600, 2.0**-500], [2.0**-600, 2.0**-574])
    emit_dot([2.0**-537, 3 * 2.0**-538], [2.0**-538, 2.0**-537])
    emit_dot([float.fromhex("0x1.8p-1022")], [0.5])
    # Overflow threshold: max + half an ulp of max ties to even, which is inf.
    emit_sum([MAX, 2.0**970])
    emit_sum([MAX, 2.0**969])
    emit_dot([2.0**600, -(2.0**600)], [2.0**600, 2.0**600])
    emit_dot([2.0**600], [2.0**600])
    # Exact zeros, signed zeros, and a vector with its own negation.
    emit_dot([0.0, -0.0], [-0.0, 5.0])
    v = [rng.uniform(-1, 1) for _ in range(40)]
    emit_sum(v + [-a for a in reversed(v)])
    # Ill-conditioned dot products (Ogita, Rump and Oishi 2005, GenDot).
    for n in (6, 39, 100):
        for _ in range(4):
            x = [rng.uniform(-1, 1) * 2.0 ** rng.randint(-60, 60) for _ in range(n)]
            y = [rng.uniform(-1, 1) * 2.0 ** rng.randint(-60, 60) for _ in range(n)]
            # Pin the last term so the exact dot nearly cancels.
            s = Fraction(0)
            for a, b in zip(x[:-1], y[:-1]):
                s += Fraction(a) * Fraction(b)
            y[-1] = 1.0
            x[-1] = -float(s) + rng.uniform(-1, 1) * 2.0**-40
            emit_dot(x, y)
    # Optimizer-shaped vectors: displacements and force differences of an
    # LJ13 relaxation (3N = 39) and a 336-atom slab (3N = 1008).
    for n, reps in ((3, 4), (39, 4), (1008, 1)):
        for _ in range(reps):
            s = [rng.gauss(0, 1e-2) for _ in range(n)]
            y = [rng.gauss(0, 1e-1) for _ in range(n)]
            emit_dot(s, y)
            emit_dot(s, s)
    # Wide dynamic range in one vector.
    x = [2.0 ** rng.randint(-1070, 1000) * rng.choice((-1, 1)) for _ in range(64)]
    emit_sum(x)
    emit_dot(x, [1.0] * len(x))


if __name__ == "__main__":
    sys.exit(main())
