import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Tactic

/-!
# Valuation Helpers
Low-level lemmas to handle p-adic valuations of integer powers.
This isolates the type-casting complexity from the main algebra.
-/

namespace PrimitiveSpectra

-- Helper 1: Valuation of p^k is k
theorem padicValRat_prime_zpow (p : ℕ) [hp : Fact p.Prime] (k : ℤ) :
    padicValRat p ((p : ℚ) ^ k) = k := by
  cases k with
  | ofNat n =>
    simp only [zpow_natCast, Int.ofNat_eq_coe]
    rw [show (p : ℚ) ^ n = ((p ^ n : ℕ) : ℚ) by norm_cast]
    rw [padicValRat.of_nat]
    rw [padicValNat.prime_pow]
  | negSucc n =>
    simp only [zpow_negSucc]
    rw [padicValRat.inv]
    rw [show (p : ℚ) ^ (n + 1) = ((p ^ (n + 1) : ℕ) : ℚ) by norm_cast]
    rw [padicValRat.of_nat]
    rw [padicValNat.prime_pow]
    rfl

-- Helper 2: Valuation of q^k is 0 if p ≠ q
theorem padicValRat_prime_other_zpow (p q : ℕ) [hp : Fact p.Prime] [hq : Fact q.Prime]
    (hneq : p ≠ q) (k : ℤ) :
    padicValRat p ((q : ℚ) ^ k) = 0 := by
  -- Pre-prove that p does not divide q
  have h_not_dvd : ¬ p ∣ q := by
    intro h_div
    -- Use basic Nat lemma: if prime p divides prime q, they are equal
    have h_eq : p = q := by
      cases (Nat.dvd_prime hq.out).mp h_div with
      | inl h1 => exact (Nat.Prime.ne_one hp.out h1).elim -- p=1 impossible
      | inr heq => exact heq -- p=q
    exact hneq h_eq

  cases k with
  | ofNat n =>
    simp only [zpow_natCast, Int.ofNat_eq_coe]
    rw [show (q : ℚ) ^ n = ((q ^ n : ℕ) : ℚ) by norm_cast]
    rw [padicValRat.of_nat]
    simp only [Nat.cast_eq_zero]
    apply padicValNat.eq_zero_of_not_dvd
    intro h_dvd_pow
    apply h_not_dvd
    exact Nat.Prime.dvd_of_dvd_pow hp.out h_dvd_pow
  | negSucc n =>
    simp only [zpow_negSucc]
    rw [padicValRat.inv]
    simp only [neg_eq_zero]
    rw [show (q : ℚ) ^ (n + 1) = ((q ^ (n + 1) : ℕ) : ℚ) by norm_cast]
    rw [padicValRat.of_nat]
    simp only [Nat.cast_eq_zero]
    apply padicValNat.eq_zero_of_not_dvd
    intro h_dvd_pow
    apply h_not_dvd
    exact Nat.Prime.dvd_of_dvd_pow hp.out h_dvd_pow

end PrimitiveSpectra
