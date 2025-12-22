import PrimitiveSpectra.Algebra.Basic
import Mathlib.Algebra.Order.Archimedean.Basic
import Mathlib.Data.Rat.Cast.Order

namespace PrimitiveSpectra
open Function

instance : Archimedean M := by
  constructor
  intro a b ha_pos
  let q_a := eval_map a
  let q_b := eval_map b
  have h_qa_gt_one : 1 < q_a := by
    rw [← eval_map_zero] at ha_pos
    -- The helper we just fixed handles the 0 case implicitly here too
    rw [eval_map]
    rw [Finsupp.prod_zero_index]
    exact ha_pos

  obtain ⟨n, h_arch⟩ := exists_pow_lt_of_lt_one h_qa_gt_one q_b
  use n
  rw [M_le, eval_map_nsmul]
  exact h_arch

end PrimitiveSpectra
