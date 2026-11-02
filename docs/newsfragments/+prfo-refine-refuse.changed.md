`prfo_refine = true` is now refused before the first force call when the client cannot refine.
That covers a client built without rgsaddle, a loaded rgsaddle library older than application binary interface (ABI) 1.15, a `prfo_refine_force` that is not positive and a `prfo_refine_max_iterations` below 1.
The check runs for a saddle search or process search with `method = min_mode` and for a nudged elastic band (NEB) job.
Before, it ran only after the search converged.
A client without rgsaddle spent 279 force calls on the dimer search of the 13-atom Lennard-Jones (LJ13) test cluster.
It then reported a potential evaluation failure.
The band of that cluster took 283 force calls before the job stopped.
Runs without `prfo_refine` are unchanged.
