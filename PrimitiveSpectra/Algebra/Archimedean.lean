import PrimitiveSpectra.Algebra.Order
import Mathlib.Algebra.Order.Archimedean.Basic
import Mathlib.Data.Rat.Cast.Order

/-!
# The Archimedean Property
...
```

The file structure is now:
```
PrimitiveSpectra/
  Algebra/
    Basic.lean      -- M, eval_map, homomorphism, positivity
    Injective.lean  -- proof of injectivity via p-adic
    Order.lean      -- LinearOrder, OrderedAddCommGroup
  Archimedean.lean  -- Archimedean instance
