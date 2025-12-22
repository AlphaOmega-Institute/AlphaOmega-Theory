import PrimitiveSpectra.Algebra.Order

/-!
# Dedekind Cuts on M
-/

namespace PrimitiveSpectra

open Set

/-- A Dedekind cut on M -/
structure Cut where
  lower : Set M
  nonempty : lower.Nonempty
  bounded : ∃ u, u ∉ lower
  downward : ∀ a ∈ lower, ∀ b, b ≤ a → b ∈ lower
  no_max : ∀ a ∈ lower, ∃ b ∈ lower, a < b

namespace Cut

/-- The upper set of a cut -/
def upper (c : Cut) : Set M := { x | x ∉ c.lower }

theorem upper_nonempty (c : Cut) : c.upper.Nonempty := c.bounded

/-- Elements of lower are less than elements of upper -/
theorem downward_not_mem (c : Cut) (a : M) (ha : a ∈ c.lower) (u : M) (hu : u ∉ c.lower) :
    a < u := by
  by_contra h
  push_neg at h
  exact hu (c.downward a ha u h)

/-- Extensionality for cuts -/
@[ext]
theorem ext {c d : Cut} (h : c.lower = d.lower) : c = d := by
  cases c; cases d; simp only [mk.injEq]; exact h

/-! ## Principal Cuts -/

/-- The principal cut at x: all elements strictly less than x -/
def principal (x : M) : Cut where
  lower := { a | a < x }
  nonempty := by
    use x - ⟨Finsupp.single ⟨2, Nat.prime_two⟩ 1⟩
    simp only [mem_setOf_eq, sub_lt_self_iff]
    rw [lt_def]
    simp only [eval_map]
    rw [Finsupp.prod_single_index]
    · simp only [zpow_one]
      exact prime_val_pos ⟨2, Nat.prime_two⟩
    · simp
  bounded := by use x; simp
  downward := by
    intro a ha b hba
    exact lt_of_le_of_lt hba ha
  no_max := by
    intro a ha
    -- Need b with a < b < x
    -- Use midpoint-like construction
    use a + ⟨Finsupp.single ⟨2, Nat.prime_two⟩ (-1)⟩
    constructor
    · simp only [mem_setOf_eq]
      calc a + ⟨Finsupp.single ⟨2, Nat.prime_two⟩ (-1)⟩
          < a + 0 := by
            rw [lt_def, eval_map_is_hom, eval_map_is_hom, eval_map_zero, mul_one]
            apply mul_lt_mul_of_pos_left _ (eval_map_pos a)
            rw [eval_map]
            rw [Finsupp.prod_single_index]
            · simp only [zpow_neg_one]
              rw [inv_lt_one_iff]
              right
              exact prime_val_pos ⟨2, Nat.prime_two⟩
            · simp
        _ = a := add_zero a
        _ < x := ha
    · rw [lt_def, eval_map_is_hom]
      calc eval_map a
          = eval_map a * 1 := (mul_one _).symm
        _ < eval_map a * eval_map ⟨Finsupp.single ⟨2, Nat.prime_two⟩ (-1)⟩⁻¹ := by
            apply mul_lt_mul_of_pos_left _ (eval_map_pos a)
            rw [eval_map, Finsupp.prod_single_index]
            · simp only [zpow_neg_one, one_lt_inv_iff]
              exact ⟨prime_val_pos _, by norm_cast; exact Nat.Prime.one_lt Nat.prime_two⟩
            · simp
        _ = eval_map a * (eval_map ⟨Finsupp.single ⟨2, Nat.prime_two⟩ (-1)⟩)⁻¹ := rfl
        _ = _ := by
            sorry -- This approach is getting complicated

theorem principal_lt {x : M} {a : M} (h : a ∈ (principal x).lower) : a < x := h

theorem principal_strictMono : StrictMono principal := by
  intro x y hxy
  sorry -- Will define order on Cut first

end Cut

end PrimitiveSpectra
