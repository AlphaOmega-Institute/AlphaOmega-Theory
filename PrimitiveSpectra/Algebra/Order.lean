import PrimitiveSpectra.Algebra.Injective
import Mathlib.Algebra.Order.Group.Defs

/-!
# Order Structure on M

The linear order and ordered group structure on M.
-/

noncomputable section

open Function

namespace PrimitiveSpectra

/-! ## Linear Order via eval_map -/

instance : LinearOrder M := LinearOrder.lift' eval_map eval_map_injective'

theorem le_def (a b : M) : a ≤ b ↔ eval_map a ≤ eval_map b := Iff.rfl

/-! ## Ordered Group Structure -/

theorem add_le_add_left_M (a b : M) (h : a ≤ b) (c : M) : c + a ≤ c + b := by
  rw [le_def] at h ⊢
  rw [eval_map_is_hom, eval_map_is_hom]
  exact mul_le_mul_of_nonneg_left h (le_of_lt (eval_map_pos c))

instance : OrderedAddCommGroup M where
  add_le_add_left := fun a b h c => add_le_add_left_M a b h c

end PrimitiveSpectra
