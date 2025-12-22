import PrimitiveSpectra.Algebra.Order
import Mathlib.Algebra.Order.Archimedean.Basic
import Mathlib.Data.Rat.Cast.Order

namespace PrimitiveSpectra

open Function

/-- For q > 1, any r is eventually exceeded by q^n -/
theorem exists_pow_gt_of_one_lt {q : ℚ} (hq : 1 < q) (r : ℚ) : ∃ n : ℕ, r < q ^ n := by
  by_cases hr : r ≤ 0
  · use 1; calc r ≤ 0 := hr
      _ < 1 := one_pos
      _ < q := hq
      _ = q ^ 1 := (pow_one q).symm
  · push_neg at hr
    have hq_pos : 0 < q := by linarith
    have hq_sub : 0 < q - 1 := by linarith
    obtain ⟨n, hn⟩ := exists_nat_gt (r / (q - 1))
    use n
    have bernoulli : ∀ m : ℕ, 1 + m * (q - 1) ≤ q ^ m := by
      intro m
      induction m with
      | zero => simp
      | succ m ih =>
        have hqm : 0 < q ^ m := pow_pos hq_pos m
        calc 1 + (m + 1) * (q - 1) = (1 + m * (q - 1)) + (q - 1) := by ring
          _ ≤ q ^ m + (q - 1) := by linarith
          _ ≤ q ^ m * q := by nlinarith
          _ = q ^ (m + 1) := (pow_succ q m).symm
    calc r < (q - 1) * n := by
          have := (div_lt_iff hq_sub).mp hn
          linarith
      _ ≤ 1 + n * (q - 1) := by linarith
      _ ≤ q ^ n := bernoulli n

instance : Archimedean M := by
  constructor
  intro a b ha_pos
  have h_one_lt : 1 < eval_map a := by
    rw [← eval_map_zero]
    rw [lt_def] at ha_pos
    exact ha_pos
  obtain ⟨n, hn⟩ := exists_pow_gt_of_one_lt h_one_lt (eval_map b)
  use n
  rw [lt_def, eval_map_nsmul]
  exact hn

end PrimitiveSpectra
