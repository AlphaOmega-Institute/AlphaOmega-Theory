import Mathlib.Algebra.Group.Basic
import Mathlib.Data.Finsupp.Basic
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.Rat.Defs
import Mathlib.Data.Rat.Lemmas
import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Tactic

/-!
# Primitive Spectra: Core Definitions
The Exponent-Vector Group M and the Evaluation Map.
-/

noncomputable section

open Function BigOperators Classical Finsupp

namespace PrimitiveSpectra

/-! ## 1. The Objects -/

def Primes : Type := { p : ℕ // p.Prime }

structure M where
  toFinsupp : Primes →₀ ℤ

namespace M

@[ext]
theorem ext {a b : M} (h : a.toFinsupp = b.toFinsupp) : a = b := by
  cases a; cases b; congr

instance : Zero M := ⟨⟨0⟩⟩
instance : Add M := ⟨fun a b => ⟨a.toFinsupp + b.toFinsupp⟩⟩
instance : Neg M := ⟨fun a => ⟨-a.toFinsupp⟩⟩
instance : Sub M := ⟨fun a b => ⟨a.toFinsupp - b.toFinsupp⟩⟩

@[simp] theorem toFinsupp_zero : (0 : M).toFinsupp = 0 := rfl
@[simp] theorem toFinsupp_add (a b : M) : (a + b).toFinsupp = a.toFinsupp + b.toFinsupp := rfl
@[simp] theorem toFinsupp_neg (a : M) : (-a).toFinsupp = -a.toFinsupp := rfl
@[simp] theorem toFinsupp_sub (a b : M) : (a - b).toFinsupp = a.toFinsupp - b.toFinsupp := rfl

instance : AddCommGroup M where
  add_assoc a b c := by ext; simp [add_assoc]
  zero_add a := by ext; simp
  add_zero a := by ext; simp
  add_comm a b := by ext; simp [add_comm]
  neg_add_cancel a := by ext; simp
  sub_eq_add_neg a b := by ext; simp [sub_eq_add_neg]
  nsmul := nsmulRec
  zsmul := zsmulRec

instance : FunLike M Primes ℤ where
  coe m := m.toFinsupp
  coe_injective' := by intro a b h; ext; exact congrFun h _

@[simp] theorem coe_toFinsupp (m : M) (p : Primes) : m.toFinsupp p = m p := rfl

end M

/-! ## 2. The Evaluation Map -/

def eval_map (v : M) : ℚ :=
  v.toFinsupp.prod (fun p exp => (p.val : ℚ) ^ exp)

/-! ## 3. Basic Properties -/

@[simp]
theorem eval_map_zero : eval_map 0 = 1 := by
  rw [eval_map, M.toFinsupp_zero]; exact Finsupp.prod_zero_index

theorem prime_val_pos (p : Primes) : (0 : ℚ) < p.val := by
  norm_cast; exact Nat.Prime.pos p.prop

theorem prime_val_ne_zero (p : Primes) : (p.val : ℚ) ≠ 0 := by
  norm_cast; exact Nat.Prime.ne_zero p.prop

theorem eval_map_is_hom (a b : M) : eval_map (a + b) = eval_map a * eval_map b := by
  rw [eval_map, eval_map, eval_map, M.toFinsupp_add]
  apply Finsupp.prod_add_index'
  · intro p; simp only [zpow_zero]
  · intro p m n; apply zpow_add₀; exact prime_val_ne_zero p

theorem eval_map_neg (v : M) : eval_map (-v) = (eval_map v)⁻¹ := by
  have h : eval_map v * eval_map (-v) = 1 := by rw [← eval_map_is_hom]; simp
  exact eq_inv_of_mul_eq_one_right h

theorem eval_map_nsmul (n : ℕ) (v : M) : eval_map (n • v) = (eval_map v) ^ n := by
  induction n with
  | zero => rw [zero_smul, pow_zero, eval_map_zero]
  | succ n ih => rw [succ_nsmul, eval_map_is_hom, ih, pow_succ]

/-! ## 4. Positivity -/

theorem zpow_prime_pos (p : Primes) (n : ℤ) : (0 : ℚ) < (p.val : ℚ) ^ n := by
  apply zpow_pos; exact prime_val_pos p

theorem eval_map_pos (v : M) : 0 < eval_map v := by
  unfold eval_map
  apply Finset.prod_pos
  intro p _
  exact zpow_prime_pos p (v.toFinsupp p)

theorem eval_map_ne_zero (v : M) : eval_map v ≠ 0 := ne_of_gt (eval_map_pos v)

end PrimitiveSpectra
