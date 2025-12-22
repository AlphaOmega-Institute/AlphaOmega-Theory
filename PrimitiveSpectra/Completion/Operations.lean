import PrimitiveSpectra.Completion.Cut

/-!
# Arithmetic Operations on Cuts
-/

namespace PrimitiveSpectra
namespace Cut

open Set

/-! ## Order on Cuts -/

instance : LE Cut where
  le c d := c.lower ⊆ d.lower

instance : LT Cut where
  lt c d := c.lower ⊂ d.lower

theorem le_def' (c d : Cut) : c ≤ d ↔ c.lower ⊆ d.lower := Iff.rfl

instance : Preorder Cut where
  le_refl c := Subset.rfl
  le_trans _ _ _ := Subset.trans
  lt_iff_le_not_le c d := by
    constructor
    · intro ⟨h, hne⟩; exact ⟨h, fun h' => hne (Subset.antisymm h h')⟩
    · intro ⟨h, hn⟩; exact ⟨h, fun heq => hn (heq ▸ Subset.rfl)⟩

instance : PartialOrder Cut where
  le_antisymm c d hcd hdc := ext (Subset.antisymm hcd hdc)

/-! ## Addition -/

def add (c d : Cut) : Cut where
  lower := { z | ∃ x ∈ c.lower, ∃ y ∈ d.lower, z = x + y }
  nonempty := by
    obtain ⟨x, hx⟩ := c.nonempty
    obtain ⟨y, hy⟩ := d.nonempty
    exact ⟨x + y, x, hx, y, hy, rfl⟩
  bounded := by
    obtain ⟨u, hu⟩ := c.bounded
    obtain ⟨v, hv⟩ := d.bounded
    use u + v
    intro ⟨x, hx, y, hy, hz⟩
    have h1 : x < u := c.downward_not_mem x hx u hu
    have h2 : y < v := d.downward_not_mem y hy v hv
    rw [hz] at h1 h2 ⊢
    exact not_lt.mpr (le_refl _) (add_lt_add h1 h2)
  downward := by
    intro z ⟨x, hx, y, hy, rfl⟩ w hwz
    use w - y + y - (x + y - w)
    sorry -- technical
  no_max := by
    intro z ⟨x, hx, y, hy, rfl⟩
    obtain ⟨x', hx', hxx'⟩ := c.no_max x hx
    exact ⟨x' + y, ⟨x', hx', y, hy, rfl⟩, add_lt_add_right hxx' y⟩

instance : Add Cut := ⟨add⟩

/-! ## Zero -/

instance : Zero Cut := ⟨principal 0⟩

/-! ## Group laws (statements) -/

theorem add_comm (c d : Cut) : c + d = d + c := by
  apply ext; ext z
  constructor <;> (rintro ⟨x, hx, y, hy, rfl⟩; exact ⟨y, hy, x, hx, add_comm x y⟩)

theorem add_assoc (a b c : Cut) : a + b + c = a + (b + c) := by
  sorry

theorem add_zero (c : Cut) : c + 0 = c := by
  sorry

theorem zero_add (c : Cut) : 0 + c = c := by
  rw [add_comm]; exact add_zero c

/-! ## Negation (hardest part) -/

def neg (c : Cut) : Cut where
  lower := { x | ∃ u ∈ c.upper, u < -x }
  nonempty := by sorry
  bounded := by sorry
  downward := by sorry
  no_max := by sorry

instance : Neg Cut := ⟨neg⟩

theorem neg_add_cancel (c : Cut) : -c + c = 0 := by
  sorry

instance : AddCommGroup Cut where
  add_assoc := add_assoc
  zero_add := zero_add
  add_zero := add_zero
  add_comm := add_comm
  neg_add_cancel := neg_add_cancel
  nsmul := nsmulRec
  zsmul := zsmulRec

end Cut
end PrimitiveSpectra
