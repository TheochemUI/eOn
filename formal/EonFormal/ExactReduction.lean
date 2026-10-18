import Mathlib

/-! # Exactness contracts behind `eonc::exact`

`client/ExactReduction.cpp` computes a dot product as an integer sum and
rounds it once. These lemmas state what that relies on, and the classical
grid lemma that makes the Ozaki splitting error-free.

* `tree_exact`: floating-point addition of numbers on a common grid
  `k * 2^q` is exact in every parenthesisation, provided the integer
  weights sum (in absolute value) to at most `2^p` and the grid sits inside
  the exponent range. No rounding mode is assumed beyond `fl x = x` for a
  representable `x`, which every IEEE mode satisfies.
* `order_independent`: two summation trees over the same integers give the
  same floating-point result. This is reproducibility.
* `slice_products_exact`: the Ozaki condition. Products of slices with
  `ba` and `bb` bits lie on the product grid; `n` of them sum exactly
  whenever `n * 2^(ba + bb) <= 2^p`.
* `product_on_grid`: a product of two significand-exponent pairs is the
  integer product on the summed exponent, the identity the accumulator
  adds.
* `rne_spec`: the rounding step of `Accumulator::round` (keep `N / 2^L`,
  add one when the half bit is set and either a sticky bit is set or the
  kept part is odd) returns the nearest integer to `N / 2^L`, ties to
  even.
* `word_capacity`: a signed 64-bit word that receives at most `2^30`
  pieces below `2^32` stays below `2^63`.
-/

namespace EonFormal

noncomputable section

/-- A binary format with `p`-bit significands and exponents in
`[emin, emax]` (subnormals are the `|m| < 2^(p-1)` values at `emin`). -/
structure Fmt where
  p : ℕ
  emin : ℤ
  emax : ℤ

/-- `x` is representable in the format. -/
def Fmt.Repr (f : Fmt) (x : ℝ) : Prop :=
  ∃ m e : ℤ, x = m * (2 : ℝ) ^ e ∧ |m| ≤ 2 ^ f.p ∧ f.emin ≤ e ∧ e ≤ f.emax

/-- A summation order: a binary tree whose leaves are integer weights on
the grid `2^q`. -/
inductive Tree
  | leaf : ℤ → Tree
  | node : Tree → Tree → Tree

namespace Tree

/-- The integer sum of the leaf weights. -/
def ints : Tree → ℤ
  | leaf k => k
  | node l r => ints l + ints r

/-- The sum of the absolute leaf weights. -/
def mass : Tree → ℤ
  | leaf k => |k|
  | node l r => mass l + mass r

/-- The number of leaves. -/
def size : Tree → ℕ
  | leaf _ => 1
  | node l r => size l + size r

/-- The leaf weights, in order. -/
def leaves : Tree → List ℤ
  | leaf k => [k]
  | node l r => leaves l ++ leaves r

/-- The exact real value of the sum. -/
def exact (q : ℤ) (t : Tree) : ℝ := (ints t : ℝ) * (2 : ℝ) ^ q

/-- The floating-point evaluation: every node rounds its sum with `fl`. -/
def eval (fl : ℝ → ℝ) (q : ℤ) : Tree → ℝ
  | leaf k => (k : ℝ) * (2 : ℝ) ^ q
  | node l r => fl (eval fl q l + eval fl q r)

theorem mass_nonneg : ∀ t : Tree, 0 ≤ mass t
  | leaf k => abs_nonneg k
  | node l r => add_nonneg (mass_nonneg l) (mass_nonneg r)

theorem abs_ints_le_mass : ∀ t : Tree, |ints t| ≤ mass t
  | leaf _ => le_refl _
  | node l r =>
    (abs_add_le _ _).trans (add_le_add (abs_ints_le_mass l) (abs_ints_le_mass r))

end Tree

open _root_.EonFormal.Tree

/-- Grid lemma: with every partial sum's integer weight bounded by `2^p`
and the grid exponent in range, rounding never acts. -/
theorem tree_exact (f : Fmt) (fl : ℝ → ℝ)
    (hfl : ∀ x, f.Repr x → fl x = x) (q : ℤ) (hq0 : f.emin ≤ q)
    (hq1 : q ≤ f.emax) :
    ∀ t : Tree, mass t ≤ 2 ^ f.p → eval fl q t = exact q t
  | leaf k, _ => by simp [eval, exact, ints]
  | node l r, h => by
    have hl : mass l ≤ 2 ^ f.p := by
      have := mass_nonneg r
      simp only [mass] at h
      linarith
    have hr : mass r ≤ 2 ^ f.p := by
      have := mass_nonneg l
      simp only [mass] at h
      linarith
    have el := tree_exact f fl hfl q hq0 hq1 l hl
    have er := tree_exact f fl hfl q hq0 hq1 r hr
    simp only [eval, el, er, exact, ints]
    have hrep : f.Repr (((ints l : ℝ)) * (2 : ℝ) ^ q + (ints r : ℝ) * (2 : ℝ) ^ q) := by
      refine ⟨ints l + ints r, q, by push_cast; ring, ?_, hq0, hq1⟩
      have := abs_ints_le_mass (node l r)
      simp only [ints, mass] at this h
      exact this.trans h
    rw [hfl _ hrep]
    push_cast
    ring

/-- Reproducibility: any two summation orders over the same integer total
agree bit for bit. -/
theorem order_independent (f : Fmt) (fl : ℝ → ℝ)
    (hfl : ∀ x, f.Repr x → fl x = x) (q : ℤ) (hq0 : f.emin ≤ q)
    (hq1 : q ≤ f.emax) (t₁ t₂ : Tree) (h₁ : mass t₁ ≤ 2 ^ f.p)
    (h₂ : mass t₂ ≤ 2 ^ f.p) (hsum : ints t₁ = ints t₂) :
    eval fl q t₁ = eval fl q t₂ := by
  rw [tree_exact f fl hfl q hq0 hq1 t₁ h₁, tree_exact f fl hfl q hq0 hq1 t₂ h₂,
    exact, exact, hsum]

theorem mass_le_size_mul (B : ℤ) :
    ∀ t : Tree, (∀ k, k ∈ leaves t → |k| ≤ B) → mass t ≤ (size t : ℤ) * B
  | leaf k, h => by simpa [mass, size, leaves] using h k (by simp [leaves])
  | node l r, h => by
    have hl := mass_le_size_mul B l (fun k hk => h k (by simp [leaves, hk]))
    have hr := mass_le_size_mul B r (fun k hk => h k (by simp [leaves, hk]))
    simp only [mass, size]
    push_cast
    linarith

/-- The Ozaki condition: `n` slice products, each a product of a `ba`-bit
and a `bb`-bit integer on the grid `2^q`, sum exactly in any order when
`n * 2^(ba + bb) <= 2^p`. -/
theorem slice_products_exact (f : Fmt) (fl : ℝ → ℝ)
    (hfl : ∀ x, f.Repr x → fl x = x) (q : ℤ) (hq0 : f.emin ≤ q)
    (hq1 : q ≤ f.emax) (ba bb : ℕ) (t : Tree)
    (hleaf : ∀ k, k ∈ leaves t →
      ∃ a b : ℤ, k = a * b ∧ |a| ≤ 2 ^ ba ∧ |b| ≤ 2 ^ bb)
    (hn : (size t : ℤ) * 2 ^ (ba + bb) ≤ 2 ^ f.p) :
    eval fl q t = exact q t := by
  apply tree_exact f fl hfl q hq0 hq1 t
  refine le_trans (mass_le_size_mul (2 ^ (ba + bb)) t ?_) hn
  intro k hk
  obtain ⟨a, b, rfl, ha, hb⟩ := hleaf k hk
  rw [abs_mul, pow_add]
  exact mul_le_mul ha hb (abs_nonneg _) (by positivity)

/-- A product of two significand-exponent pairs is the integer product of
the significands on the summed exponent, with no rounding. -/
theorem product_on_grid (mx my : ℤ) (ex ey : ℤ) :
    ((mx : ℝ) * (2 : ℝ) ^ ex) * ((my : ℝ) * (2 : ℝ) ^ ey) =
      ((mx * my : ℤ) : ℝ) * (2 : ℝ) ^ (ex + ey) := by
  rw [zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]
  push_cast
  ring

/-- Two 53-bit significands multiply below `2^106`: the product fits the
five 32-bit accumulator pieces after any shift below 32. -/
theorem product_width (mx my : ℕ) (hx : mx < 2 ^ 53) (hy : my < 2 ^ 53)
    (s : ℕ) (hs : s < 32) : mx * my * 2 ^ s < 2 ^ (32 * 5) := by
  have h1 : mx * my < 2 ^ 106 := by
    calc mx * my < 2 ^ 53 * 2 ^ 53 :=
          Nat.mul_lt_mul'' hx hy
      _ = 2 ^ 106 := by norm_num
  have h2 : 2 ^ s ≤ 2 ^ 31 := Nat.pow_le_pow_right (by norm_num) (by omega)
  calc mx * my * 2 ^ s < 2 ^ 106 * 2 ^ 31 := by
        exact Nat.mul_lt_mul_of_lt_of_le h1 h2 (by positivity)
    _ ≤ 2 ^ (32 * 5) := by norm_num

/-- The rounding step of `Accumulator::round` on the integer `N` with `L`
dropped bits, written with `D = 2^L = 2 H`. -/
def rne (N H : ℕ) : ℕ :=
  let D := 2 * H
  let q := N / D
  let r := N % D
  if r / H = 1 ∧ (r % H ≠ 0 ∨ q % 2 = 1) then q + 1 else q

/-- `rne` is round to nearest, ties to even: the distance to `N / D` is at
most one half, and exactly one half only at an even result. Stated
without division as `|2 (N - z D)| <= D`. -/
theorem rne_spec (N H : ℕ) (hH : 0 < H) :
    let D := 2 * H
    let z := rne N H
    |2 * ((N : ℤ) - z * D)| ≤ D ∧
      (|2 * ((N : ℤ) - z * D)| = D → z % 2 = 0) := by
  intro D z
  have hD : 0 < D := by omega
  have hdiv := Nat.div_add_mod N D
  have hr : N % D < D := Nat.mod_lt N hD
  have hhalf := Nat.div_add_mod (N % D) H
  have hs : N % D % H < H := Nat.mod_lt _ hH
  -- the half bit is 0 or 1 because r < 2 H
  have hb : N % D / H ≤ 1 := by
    have : N % D / H < 2 := (Nat.div_lt_iff_lt_mul hH).2 (by omega)
    omega
  set q := N / D with hq
  set r := N % D with hrdef
  set h := r / H with hh
  set s := r % H with hsdef
  have hN : (N : ℤ) = q * D + r := by
    have : (N : ℤ) = ((D * q + r : ℕ) : ℤ) := by rw [hdiv]
    rw [this]; push_cast; ring
  have hR : (r : ℤ) = h * H + s := by
    have : (r : ℤ) = ((H * h + s : ℕ) : ℤ) := by rw [hhalf]
    rw [this]; push_cast; ring
  have hDz : (D : ℤ) = 2 * H := by simp [D]
  simp only [z, rne]
  split_ifs with hc
  · -- rounded up: N - (q + 1) D = r - D
    have e : (N : ℤ) - ((q + 1 : ℕ) : ℤ) * D = r - D := by
      rw [hN]; push_cast; ring
    rw [e]
    obtain ⟨h1, hst⟩ := hc
    have hh1 : (h : ℤ) = 1 := by exact_mod_cast h1
    rcases hst with hst | hst
    · have hs0 : (0 : ℤ) < s := by exact_mod_cast Nat.pos_of_ne_zero hst
      constructor
      · rw [abs_le]; constructor <;> nlinarith
      · intro habs
        exfalso
        rw [abs_eq (by positivity)] at habs
        rcases habs with habs | habs <;> nlinarith
    · constructor
      · rw [abs_le]; constructor <;> nlinarith
      · intro _
        omega
  · -- kept: N - q D = r
    have e : (N : ℤ) - (q : ℤ) * D = r := by rw [hN]; ring
    rw [e]
    have hr0 : (0 : ℤ) ≤ r := by positivity
    rw [abs_of_nonneg (by positivity)]
    push Not at hc
    rcases Nat.le_one_iff_eq_zero_or_eq_one.1 hb with h0 | h1
    · have hh0 : (h : ℤ) = 0 := by exact_mod_cast h0
      constructor
      · nlinarith
      · intro habs; exfalso; nlinarith
    · obtain ⟨hs0, hq0⟩ := hc h1
      have hh1 : (h : ℤ) = 1 := by exact_mod_cast h1
      have hs0' : (s : ℤ) = 0 := by exact_mod_cast hs0
      constructor
      · nlinarith
      · intro _; omega

theorem list_sum_bound :
    ∀ pieces : List ℤ, (∀ x ∈ pieces, |x| < 2 ^ 32) →
      |pieces.sum| ≤ pieces.length * 2 ^ 32
  | [], _ => by simp
  | x :: xs, hp => by
    have hx := hp x (by simp)
    have ih := list_sum_bound xs (fun y hy => hp y (by simp [hy]))
    simp only [List.sum_cons, List.length_cons]
    push_cast
    calc |x + xs.sum| ≤ |x| + |xs.sum| := abs_add_le _ _
      _ ≤ 2 ^ 32 + xs.length * 2 ^ 32 := by linarith
      _ = (xs.length + 1) * 2 ^ 32 := by ring

/-- A word that starts below `2^32` and receives at most `2^30` signed
pieces below `2^32` stays inside a signed 64-bit integer. -/
theorem word_capacity (w : ℤ) (pieces : List ℤ) (hw : |w| < 2 ^ 32)
    (hp : ∀ x ∈ pieces, |x| < 2 ^ 32) (hn : pieces.length ≤ 2 ^ 30) :
    |w + pieces.sum| < 2 ^ 63 := by
  have hsum := list_sum_bound pieces hp
  have hlen : (pieces.length : ℤ) * 2 ^ 32 ≤ 2 ^ 30 * 2 ^ 32 := by
    have : (pieces.length : ℤ) ≤ 2 ^ 30 := by exact_mod_cast hn
    nlinarith
  calc |w + pieces.sum| ≤ |w| + |pieces.sum| := abs_add_le _ _
    _ < 2 ^ 32 + 2 ^ 30 * 2 ^ 32 := by linarith
    _ ≤ 2 ^ 63 := by norm_num

end

end EonFormal
