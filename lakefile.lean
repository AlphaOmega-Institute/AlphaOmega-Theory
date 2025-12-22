import Lake
open Lake DSL

package «PrimitiveSpectra» {
  moreLeanArgs := #["-DautoImplicit=false"]
}

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.15.0"

@[default_target]
lean_lib «PrimitiveSpectra» {
}
