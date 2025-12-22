import PrimitiveSpectra.Algebra.Basic
import PrimitiveSpectra.Algebra.Valuation -- Import the clean helpers
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Tactic

namespace PrimitiveSpectra

open Finsupp BigOperators

instance : DecidableEq Primes := Subtype.instDecidableEq

/-! ### 1. Helper Lemma: Valuation of a Product -/

theorem padicValRat_finset_prod {α : Type*} [DecidableEq α] (s : Finset α) (f : α → ℚ)
    (p : ℕ) [Fact p.Prime] (h_nonzero : ∀ x ∈ s, f x ≠ 0) :
    padicValRat p (∏ x in s, f x) = ∑ x in s, padicValRat p (f x) := by
  induction s using Finset.induction_on with
  | empty => simp only [Finset.prod_empty, Finset.sum_empty, padicValRat.one]
  | insert h_not_mem ih =>
    rw [Finset.prod_insert h_not_mem, Finset.sum_insert h_not_mem]
    rw [padicValRat.mul]
    · rw [ih]
      · intro x hx; exact h_nonzero x (Finset.mem_insert_of_mem hx)
    · exact h_nonzero _ (Finset.mem_insert_self _ _)
    · apply Finset.prod_ne_zero_iff.mpr
      intro x hx
      exact h_nonzero x (Finset.mem_insert_of_mem hx)

/-! ### 2. The Theorem -/

theorem padicValRat_eval_map (v : M) (p : Primes) :
    padicValRat p.val (eval_map v) = v p := by
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  rw [eval_map, Finsupp.prod, padicValRat_finset_prod]

  by_cases hp : p ∈ v.toFinsupp.support
  · rw [Finset.sum_eq_single p]
    · simp only [M.coe_toFinsupp]
      -- CLEAN: We use our specific helper from Valuation.lean
      rw [padicValRat_prime_zpow]
    · intro q hq hneq
      simp only [M.coe_toFinsupp]
      -- CLEAN: We use our specific helper
      have h_ne_val : p.val ≠ q.val := by
        intro h; apply hneq; apply Subtype.ext; exact h.symm

      haveI : Fact (Nat.Prime q.val) := ⟨q.prop⟩
      rw [padicValRat_prime_other_zpow p.val q.val h_ne_val]
    · intro h; contradiction
  · rw [Finset.sum_eq_zero]
    · simp at hp; rw [M.coe_toFinsupp, Finsupp.not_mem_support_iff.mp hp]
    · intro q hq
      simp only [M.coe_toFinsupp]
      have h_ne_val : p.val ≠ q.val := by
        intro h; apply hp; rw [←h] at hq; exact hq -- Logic: if p=q, p is in support

      haveI : Fact (Nat.Prime q.val) := ⟨q.prop⟩
      rw [padicValRat_prime_other_zpow p.val q.val h_ne_val]

  -- Nonzero requirement
  · intro q hq
    apply zpow_ne_zero
    norm_cast; exact Nat.Prime.ne_zero q.prop

theorem eval_map_injective' : Function.Injective eval_map := by
  intro a b h
  apply M.ext
  ext p
  dsimp [M.coe_toFinsupp]
  haveI : Fact (Nat.Prime p.val) := ⟨p.prop⟩
  rw [← padicValRat_eval_map a p, h, padicValRat_eval_map b p]

end PrimitiveSpectra
