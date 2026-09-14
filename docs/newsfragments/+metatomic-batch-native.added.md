Metatomic `forceBatch` tries a single `model.forward` over N systems and falls back to sequential `force()` if that path throws.
