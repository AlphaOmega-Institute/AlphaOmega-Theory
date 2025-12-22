import Mathlib.Algebra.Group.Basic
import Mathlib.Data.Finsupp.Basic
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.Rat.Defs
import Mathlib.Data.Rat.Lemmas
import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Tactic

/-!
# Primitive Spectra: The Arithmetic Substrate
Definitions of the Exponent-Vector Group M and the Evaluation Map.
-/

noncomputable section

open Function BigOperators

namespace PrimitiveSpectra

-- 1. THE OBJECTS
def Primes : Type := { p : ℕ // p.Prime }
def M : Type := Primes →₀ ℤ
instance : AddCommGroup M := Finsupp.instAddCommGroup

-- 2. THE BRIDGE
def eval_map (v : M) : ℚ :=
  v.prod (fun p exp => (p.val : ℚ) ^ exp)

-- 3. THE STRUCTURAL ISOMORPHISM
theorem eval_map_is_hom (a b : M) : eval_map (a + b) = eval_map a * eval_map b := by
  rw [eval_map, eval_map, eval_map]
  apply Finsupp.prod_add_index'
  · intro p
    simp only [zpow_zero]
  · intro p m n
    apply zpow_add₀
    norm_cast
    exact (Nat.Prime.ne_zero p.prop)

-- 4. THE ORDERING
def M_le (a b : M) : Prop := eval_map a ≤ eval_map b
instance : LE M := ⟨M_le⟩

-- 5. HELPER LEMMA (Final Version)
theorem eval_map_nsmul (n : ℕ) (v : M) : eval_map (n • v) = (eval_map v) ^ n := by
  induction' n with n ih
  · -- Case 0: 0 • v = 0, and anything^0 = 1.
    rw [zero_smul, pow_zero]
    -- Prove eval_map 0 = 1 explicitly using the Finsupp library
    rw [eval_map]
    exact Finsupp.prod_zero_index
  · -- Case n + 1
    rw [succ_nsmul, eval_map_is_hom, ih]
    ring

end PrimitiveSpectra
