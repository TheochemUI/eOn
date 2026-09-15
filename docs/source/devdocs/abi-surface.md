# Public C++ surface

The installed `include/eon` tree still ships internal headers. The
Job factory has landed. `Parameters` now has a private load-state
`Impl` (`last_load_source` / `last_load_error`); option-group layout
stays in the installed header and is **not** ABI-stable. `Matter`
still exposes Eigen members in the header (accessors return Eigen
types). Treat only these as the extension contract:

| Type | Pure virtuals | Notes |
|---|---|---|
| `Potential` | `force(...)` | Raw C arrays stay the Fortran/FFI adapter. `get_ef` is Eigen. |
| `Job` | `run()` | `runFromMatter` is extra, not required. |
| `Optimizer` | `step`, `run` | New methods must be defaulted virtuals, not new pures. |

Capability queries (`isThreadSafe`, `layoutFlags`, `supportsBatchEvaluation`)
are defaulted virtuals. Do not add another `= 0` without a towncrier
breaking fragment.

CI should keep grepping that those three pures still exist and that new
pures are not added silently.
