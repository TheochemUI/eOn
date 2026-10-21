# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy>=1.12"]
# ///
"""Derivations behind client/ExactReduction.cpp, checked with sympy.

Run with:  uv run --with sympy client/validation/exact_reduction_derivations.py

1. The Ozaki splitting condition (Ozaki, Ogita, Oishi and Rump 2012).
   With p-bit significands, a slice of x is a multiple of 2^(tau+beta-p)
   no larger than 2^tau: an integer weight of at most p - beta bits. A
   product of two slices then carries at most 2 (p - beta) bits on the
   product grid, and n of them sum without rounding iff
   n 2^(2 (p - beta)) <= 2^p, i.e. beta >= (p + log2 n) / 2.
2. The portable 64 x 64 -> 128 bit product: the four 32-bit partial
   products and the middle carry recombine to a * b, and the middle
   column stays below 3 * 2^32, so it fits in 64 bits.
3. The five 32-bit pieces of (hi, lo) * 2^s that Accumulator::add writes
   sum back to the shifted product.
4. Why Eigen 3.4 and 5.0 give different bits for the same dot product of
   length 39 under SSE2 (two doubles per packet): the two releases build
   different summation trees, so the first-order rounding errors carry
   different coefficients. The exact path has none.
"""

from __future__ import annotations

import math

import sympy as sp


def ozaki_condition() -> None:
    p, beta, n = sp.symbols("p beta n", positive=True)
    # n products of (p - beta)-bit slices on one grid: the integer weight
    # of the sum is at most n 2^(2 (p - beta)); exact iff <= 2^p.
    log_lhs = sp.expand_log(sp.log(n * 2 ** (2 * (p - beta)), 2), force=True)
    assert sp.simplify(log_lhs - (sp.log(n, 2) + 2 * (p - beta))) == 0
    beta_min = sp.solve(sp.Eq(log_lhs, p), beta)[0]
    expect = (p + sp.log(n, 2)) / 2
    assert sp.simplify(beta_min - expect) == 0, beta_min
    print(f"1. beta_min = {sp.simplify(expect)}")
    for nn in (3, 39, 1008, 30000):
        b = math.ceil((53 + math.log2(nn)) / 2)
        bits = 53 - b
        assert nn * 2 ** (2 * bits) <= 2**53
        print(f"   n = {nn:6d}: beta = {b}, {bits} bits per slice")


def portable_mul64() -> None:
    a0, a1, b0, b1 = sp.symbols("a0 a1 b0 b1", integer=True, nonnegative=True)
    W = sp.Integer(2) ** 32
    p00, p01, p10, p11 = a0 * b0, a0 * b1, a1 * b0, a1 * b1
    # Column identity of the schoolbook product.
    assert sp.expand((a1 * W + a0) * (b1 * W + b0) - (p00 + (p01 + p10) * W + p11 * W**2)) == 0
    # Digits: p00 = h00 W + l00, p01 = h01 W + l01, p10 = h10 W + l10,
    # mid = h00 + l01 + l10 = hm W + lm.
    h00, l00, h01, l01, h10, l10, hm = sp.symbols(
        "h00 l00 h01 l01 h10 l10 hm", integer=True, nonnegative=True
    )
    q11 = sp.Symbol("p11", integer=True, nonnegative=True)
    mid = h00 + l01 + l10
    lm = mid - hm * W
    lo = l00 + lm * W  # (p00 & m) | (mid << 32), only mid's low digit kept
    hi = q11 + h01 + h10 + hm
    columns = (h00 * W + l00) + ((h01 * W + l01) + (h10 * W + l10)) * W + q11 * W**2
    assert sp.expand(hi * W**2 + lo - columns) == 0
    # Each digit is below 2^32, so mid <= 3 (2^32 - 1) fits in 64 bits.
    assert 3 * (2**32 - 1) < 2**64
    print("2. hi 2^64 + lo = a b; middle column < 3 * 2^32")


def accumulator_pieces() -> None:
    W = sp.Integer(2) ** 32
    d = sp.symbols("d0:6", integer=True, nonnegative=True)
    # (hi, lo) << s as 32-bit digits d0..d4 (d5 is zero because the
    # shifted product is below 2^138 < 2^160).
    shifted = sum(d[k] * W**k for k in range(5))
    t0 = d[0] + d[1] * W
    t1 = d[2] + d[3] * W
    t2 = d[4]
    pieces = [t0 - d[1] * W, (t0 - d[0]) / W, t1 - d[3] * W, (t1 - d[2]) / W, t2]
    recombined = sum(pieces[k] * W**k for k in range(5))
    assert sp.simplify(recombined - shifted) == 0
    assert 2**53 * 2**53 * 2**31 < 2**160
    print("3. five pieces of (hi, lo) << s recombine; 2^106 * 2^31 < 2^160")


def eigen_trees(n: int = 39, packet: int = 2) -> None:
    # Eigen 3.4 redux (Redux.h): two packet accumulators over pairs of
    # packets, one trailing packet, horizontal add, then scalars.
    def eigen34() -> list[list[int]]:
        aligned2 = (n // (2 * packet)) * (2 * packet)
        aligned = (n // packet) * packet
        acc0 = [list(range(0, packet))]
        acc1 = [list(range(packet, 2 * packet))]
        for i in range(2 * packet, aligned2, 2 * packet):
            acc0.append(list(range(i, i + packet)))
            acc1.append(list(range(i + packet, i + 2 * packet)))
        groups = [sum(acc0, []) + sum(acc1, [])]
        if aligned > aligned2:
            groups[0] += list(range(aligned2, aligned))
        return groups + [[i] for i in range(aligned, n)]

    # Eigen 5.0 inner_product_impl (InnerProduct.h): four accumulators.
    def eigen50() -> list[list[int]]:
        quad = (n // (4 * packet)) * (4 * packet)
        pend = (n // packet) * packet
        accs: list[list[int]] = [[] for _ in range(4)]
        for i in range(0, quad, 4 * packet):
            for a in range(4):
                accs[a] += list(range(i + a * packet, i + (a + 1) * packet))
        for r, i in enumerate(range(quad, pend, packet)):
            accs[r] += list(range(i, i + packet))
        return [sum(accs, [])] + [[i] for i in range(pend, n)]

    # Rounding errors at first order: each addition into an accumulator
    # lane contributes delta * (running lane sum). Count, per term, how
    # many additions it passes through before the final scalar.
    def depth(order: str) -> list[int]:
        if order == "3.4":
            lanes = {}
            aligned2 = (n // (2 * packet)) * (2 * packet)
            aligned = (n // packet) * packet
            for i in range(aligned):
                acc = (i // packet) % 2 if i < aligned2 else 0
                lanes.setdefault((acc, i % packet), []).append(i)
        else:
            lanes = {}
            quad = (n // (4 * packet)) * (4 * packet)
            pend = (n // packet) * packet
            for i in range(pend):
                acc = (i // packet) % 4 if i < quad else (i - quad) // packet
                lanes.setdefault((acc, i % packet), []).append(i)
        return sorted(len(v) for v in lanes.values())

    e34, e50 = eigen34(), eigen50()
    assert sorted(sum(e34, [])) == sorted(sum(e50, [])) == list(range(n))
    d34, d50 = depth("3.4"), depth("5.0")
    assert d34 != d50
    print(
        f"4. n = {n}, packet = {packet}: lane lengths Eigen 3.4 {d34}, "
        f"Eigen 5.0 {d50}: different trees, different roundings"
    )


if __name__ == "__main__":
    ozaki_condition()
    portable_mul64()
    accumulator_pieces()
    eigen_trees()
