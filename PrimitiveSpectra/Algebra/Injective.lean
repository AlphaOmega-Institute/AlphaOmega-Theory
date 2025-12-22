import PrimitiveSpectra.Algebra.Basic
import Mathlib.NumberTheory.Padics.PadicVal.Basic

/-!
# Injectivity of eval_map

Proof that eval_map is injective using p-adic valuations.
This encodes the Fundamental Theorem of Arithmetic.
-/

namespace PrimitiveSpectra

open Finsupp

/-! ## p-adic valuation lemmas -/

/-- p-adic valuation of p itself -/
theorem padicValRat_self (p : Primes) :
    padicValRat p.val (p.val : ℚ) = 1 := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  exact padicValRat.self p.prop.one_lt

/-- p-adic valuation of p^n for natural n -/
theorem padicValRat_self_pow (p : Primes) (n : ℕ) :
    padicValRat p.val ((p.val : ℚ) ^ n) = n := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  induction n with
  | zero => simp [padicValRat.one]
  | succ n ih =>
    rw [pow_succ, mul_comm]
    rw [padicValRat.mul (prime_val_ne_zero p) (pow_ne_zero n (prime_val_ne_zero p))]
    rw [ih, padicValRat_self p]
    simp only [Nat.cast_succ]
    ring

/-- p-adic valuation of p^n for integer n -/
theorem padicValRat_self_zpow (p : Primes) (n : ℤ) :
    padicValRat p.val ((p.val : ℚ) ^ n) = n := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  cases n with
  | ofNat m =>
    simp only [zpow_natCast, Int.ofNat_eq_coe]
    exact_mod_cast padicValRat_self_pow p m
  | negSucc m =>
    simp only [zpow_negSucc]
    have hne : (p.val : ℚ) ^ (m + 1) ≠ 0 := pow_ne_zero _ (prime_val_ne_zero p)
    rw [padicValRat.inv hne]
    simp only [neg_inj]
    exact_mod_cast padicValRat_self_pow p (m + 1)

/-- If p and q are distinct primes, p does not divide q -/
theorem prime_not_dvd_prime (p q : Primes) (hne : p ≠ q) : ¬(p.val ∣ q.val) := by
  intro hdvd
  have h := Nat.Prime.eq_one_or_self_of_prime_of_dvd q.prop p.val p.prop hdvd
  rcases h with rfl | rfl
  · exact Nat.Prime.ne_one p.prop rfl
  · exact hne rfl

/-- p-adic valuation of a different prime q is zero -/
theorem padicValRat_other (p q : Primes) (hne : p ≠ q) :
    padicValRat p.val (q.val : ℚ) = 0 := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  rw [padicValRat.of_nat]
  simp only [Int.natCast_eq_zero]
  exact padicValNat.eq_zero_of_not_dvd (prime_not_dvd_prime p q hne)

/-- p-adic valuation of q^n for natural n is zero when q ≠ p -/
theorem padicValRat_other_pow (p q : Primes) (hne : p ≠ q) (n : ℕ) :
    padicValRat p.val ((q.val : ℚ) ^ n) = 0 := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  induction n with
  | zero => simp [padicValRat.one]
  | succ n ih =>
    rw [pow_succ, mul_comm]
    rw [padicValRat.mul (prime_val_ne_zero q) (pow_ne_zero n (prime_val_ne_zero q))]
    rw [ih, padicValRat_other p q hne]
    ring

/-- p-adic valuation of q^n is zero when q ≠ p -/
theorem padicValRat_other_zpow (p q : Primes) (hne : p ≠ q) (n : ℤ) :
    padicValRat p.val ((q.val : ℚ) ^ n) = 0 := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  cases n with
  | ofNat m =>
    simp only [zpow_natCast, Int.ofNat_eq_coe]
    exact_mod_cast padicValRat_other_pow p q hne m
  | negSucc m =>
    simp only [zpow_negSucc]
    have hne' : (q.val : ℚ) ^ (m + 1) ≠ 0 := pow_ne_zero _ (prime_val_ne_zero q)
    rw [padicValRat.inv hne']
    simp only [neg_eq_zero]
    exact_mod_cast padicValRat_other_pow p q hne (m + 1)

/-! ## Product lemma -/

/-- p-adic valuation is additive over products -/
theorem padicValRat_prod {ι : Type*} (s : Finset ι) (f : ι → ℚ) (hf : ∀ i ∈ s, f i ≠ 0)
    (p : ℕ) [Fact p.Prime] :
    padicValRat p (∏ i in s, f i) = ∑ i in s, padicValRat p (f i) := by
  induction s using Finset.induction_on with
  | empty => simp [padicValRat.one]
  | insert ha ih =>
    rw [Finset.prod_insert ha, Finset.sum_insert ha]
    have h1 : f _ ≠ 0 := hf _ (Finset.mem_insert_self _ _)
    have h2 : ∀ i ∈ _, f i ≠ 0 := fun i hi => hf i (Finset.mem_insert_of_mem hi)
    have h3 : ∏ i in _, f i ≠ 0 := Finset.prod_ne_zero_iff.mpr h2
    rw [padicValRat.mul h1 h3, ih h2]

/-! ## Main theorem -/

/-- The p-adic valuation extracts the exponent at prime p -/
theorem padicValRat_eval_map (v : M) (p : Primes) :
    padicValRat p.val (eval_map v) = v p := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  unfold eval_map
  rw [Finsupp.prod]
  rw [padicValRat_prod _ _ (fun q _ => ne_of_gt (zpow_prime_pos q (v.toFinsupp q)))]
  by_cases hp : p ∈ v.toFinsupp.support
  · rw [← Finset.add_sum_erase _ _ hp]
    simp only [M.coe_toFinsupp, padicValRat_self_zpow]
    suffices h : ∑ x ∈ v.toFinsupp.support.erase p, padicValRat (↑p) (↑↑x ^ v.toFinsupp x) = 0 by
      rw [h]; ring
    apply Finset.sum_eq_zero
    intro q hq
    have hne : q ≠ p := Finset.ne_of_mem_erase hq
    exact padicValRat_other_zpow p q hne.symm (v.toFinsupp q)
  · have hvp : v p = 0 := Finsupp.not_mem_support_iff.mp hp
    rw [hvp]
    apply Finset.sum_eq_zero
    intro q hq
    have hne : q ≠ p := fun h => hp (h ▸ hq)
    exact padicValRat_other_zpow p q hne.symm (v.toFinsupp q)

/-- Injectivity of eval_map - the Fundamental Theorem of Arithmetic -/
theorem eval_map_injective' : Function.Injective eval_map := by
  intro a b hab
  ext p
  have ha := padicValRat_eval_map a p
  have hb := padicValRat_eval_map b p
  rw [hab] at ha
  exact ha.symm.trans hb

end PrimitiveSpectra
