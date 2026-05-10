The Python orchestrator now reads termination codes through
``eon.status_codes.SaddleStatus`` / ``MinimizationStatus``
(``IntEnum`` mirrors of ``client/StatusTypes.h``) instead of bare
integer literals. The wire format is unchanged -- ``results.dat``
still carries integer ``termination_reason`` values -- but
``akmcstate``, ``explorer``, and ``basinhopping`` now compare against
named members and route bad-saddle labels through a single source of
truth. Two latent bugs surface and are fixed: ``register_bad_saddle``
no longer ``IndexError``s on codes 20-22 (DimerLostMode,
DimerRestoredBest, BadArtnError), and the ``"nonlocal abort"`` slug
typo collapses to ``"nonlocal_abort"``.
