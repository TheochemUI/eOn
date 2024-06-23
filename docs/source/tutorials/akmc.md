# AKMC Tutorial

You can find a sample setup for the AKMC method in the examples/akmc-pt
directory of your copy of eon. The system is a Pt heptamer island on a Pt(111)
surface:

```{figure} ../fig/akmc-1.png
---
alt: Pt heptamer island on Pt(111)
class: full-width
align: center
---
The Pt heptamer island on a Pt(111) surface.
```


After following the [installation instructions](project:../install/index.md):

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit 
➜ python -m eon.server
Eon version 1d9860cf
State list path does not exist; Creating: .//states/
Registering results
Processed 0 results
Queue contains 0 searches
Making 2 process searches
Job finished: .//jobs/scratch/0_0
Job finished: .//jobs/scratch/0_1
Created 2 searches
Currently in state 0 with confidence 0.000000
```

Each time you run EON, the server registers results from the previous execution.
The first time the EON server is run, there are no results to register. Two
saddle searches are executed, but their results will only be registered on the
next execution of EON. You can see this by looking at the `search_results.txt`
file in your states directory.

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit 
➜ cat states/0/search_results.txt
    wuid       type    barrier   max-dist    sad-fcs   mins-fcs   pref-fcs    result
-------------------------------------------------------------------------------------
```

It's empty. We can register the results of the first execution of EON by running
the server again

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit 
➜ python -m eon.server
Eon version 1d9860cf
Registering results
Found new lowest barrier 0.987074 for state 0 (type: random)
Found new barrier 0.987074 for state 0 (type: random)
Found new barrier 1.196473 for state 0 (type: random)
Processed 2 results
Queue contains 0 searches
Making 2 process searches
Job finished: .//jobs/scratch/0_2
Job finished: .//jobs/scratch/0_3
Created 2 searches
Currently in state 0 with confidence 0.000000
```

Now some results should be present in the table (these will vary):

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit took 10s 
➜ cat states/0/search_results.txt
    wuid       type    barrier   max-dist    sad-fcs   mins-fcs   pref-fcs    result
-------------------------------------------------------------------------------------
       0     random    0.98707    0.00000        282        339       1044    good-0
       1     random    1.19647    0.00000        630        375       1287    good-1
```

The `search_results.txt` file provides information about how a particular job
went. The data collected is contained in the ``processtable`` file in the state
directory:

```{code-block} bash

EONgit/examples/akmc-pt on main via 🅒 eongit 
➜ cat states/0/processtable 
 proc #    saddle energy   prefactor   product   product energy product prefactor  barrier         rate repeats
      0      -1774.80409 7.32913e+12        -1      -1775.00608       2.64735e+12  0.98707  1.91855e-04       0
      1      -1774.59469 1.48200e+13        -1      -1775.09375       4.55267e+12  1.19647  1.17764e-07       0
```

You can find the files relevant to a process in the `procdata` directory of a
given state:

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit 
➜ ls states/0/procdata 
mode_0.dat  product_0.con  reactant_0.con  results_0.dat  saddle_0.con  stdout_0.dat
mode_1.dat  product_1.con  reactant_1.con  results_1.dat  saddle_1.con  stdout_1.dat
```

The data relevant to process 0 is in the process table. The figures below show
the reactant, saddle, and product configurations for this process:

```{figure} ../fig/akmc-2.png
---
alt: The reaction path critical points
class: full-width
align: center
---
Left to right, reactant, saddle and product configurations.
```

Subsequent runs of EON will show an increasing confidence as the event table (in
the specified energy window) becomes complete:

```{code-block} bash
currently in state 0 with confidence 0.313698
...
currently in state 0 with confidence 0.641981
...
currently in state 0 with confidence 0.675757
```

We need to explore a larger part of the phase space until a KMC step can be
taken. To iterate over states quickly, consider a `for` loop in `bash`:

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit took 7s 
➜ for i in {0..60}; do python -m eon.server; done
```

Eventually, the simulation will reach the required confidence and take a KMC
step to the next state:

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit 
➜ python -m eon.server
Registering results
Processed 4 results
Cancelled 0 workunits from state 0
register_process: new product state
register_process: into eq rate test
Meantime for Step 1:  0.0003931769704500849
KMC step from state 0 through process 2 to state 1 
currently in state 1 with confidence 0.000000
```

This will then be reflected in the `dynamics.txt` file of the simulation directory:

```{code-block} bash
EONgit/examples/akmc-pt on main via 🅒 eongit 
➜ cat dynamics.txt
 step-number   reactant-id    process-id    product-id     step-time    total-time       barrier          rate        energy
-----------------------------------------------------------------------------------------------------------------------------
           0             0             2             1  1.531416e-03  1.531416e-03      0.601060  9.861450e+02  -1775.791160
```
