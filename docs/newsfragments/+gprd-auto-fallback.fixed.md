`with_gprd=auto` leaves the GP dimer off when the `gpr_optim` fetch fails. `with_gprd=enabled` still stops configuration.

A build directory that stored `with_gprd` as `true` or `false` rejects `meson setup --reconfigure` (`Option "with_gprd" value auto is not boolean`). `python scripts/migrate_with_gprd_option.py <builddir>` maps `true` to `enabled` and `false` to `disabled`, and the reconfigure then runs.
