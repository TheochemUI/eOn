# {py:mod}`eon.akmcstate`

```{py:module} eon.akmcstate
```

```{autodoc2-docstring} eon.akmcstate
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`AKMCState <eon.akmcstate.AKMCState>`
  - ```{autodoc2-docstring} eon.akmcstate.AKMCState
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`lambertw <eon.akmcstate.lambertw>`
  - ```{autodoc2-docstring} eon.akmcstate.lambertw
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.akmcstate.logger>`
  - ```{autodoc2-docstring} eon.akmcstate.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.akmcstate.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.akmcstate.logger
```

````

`````{py:class} AKMCState(statepath, statenumber, statelist, previous_state_num=-1, reactant_path=None, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.akmcstate.AKMCState

Bases: {py:obj}`eon.state.State`

```{autodoc2-docstring} eon.akmcstate.AKMCState
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.akmcstate.AKMCState.__init__
```

````{py:attribute} processtable_head_fmt
:canonical: eon.akmcstate.AKMCState.processtable_head_fmt
:value: >
   '%7s %16s %11s %9s %16s %17s %8s %12s %7s\n'

```{autodoc2-docstring} eon.akmcstate.AKMCState.processtable_head_fmt
```

````

````{py:attribute} processtable_header
:canonical: eon.akmcstate.AKMCState.processtable_header
:value: >
   None

```{autodoc2-docstring} eon.akmcstate.AKMCState.processtable_header
```

````

````{py:attribute} processtable_line
:canonical: eon.akmcstate.AKMCState.processtable_line
:value: >
   '%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n'

```{autodoc2-docstring} eon.akmcstate.AKMCState.processtable_line
```

````

````{py:method} find_repeat(saddle_file, barrier)
:canonical: eon.akmcstate.AKMCState.find_repeat

```{autodoc2-docstring} eon.akmcstate.AKMCState.find_repeat
```

````

````{py:method} add_process(result, superbasin=None)
:canonical: eon.akmcstate.AKMCState.add_process

```{autodoc2-docstring} eon.akmcstate.AKMCState.add_process
```

````

````{py:method} append_search_result(result, comment, superbasin)
:canonical: eon.akmcstate.AKMCState.append_search_result

```{autodoc2-docstring} eon.akmcstate.AKMCState.append_search_result
```

````

````{py:method} get_ratetable(superbasin=None)
:canonical: eon.akmcstate.AKMCState.get_ratetable

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_ratetable
```

````

````{py:method} get_process_table()
:canonical: eon.akmcstate.AKMCState.get_process_table

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_process_table
```

````

````{py:method} get_relevant_procids(superbasin=None)
:canonical: eon.akmcstate.AKMCState.get_relevant_procids

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_relevant_procids
```

````

````{py:method} get_confidence(superbasin=None)
:canonical: eon.akmcstate.AKMCState.get_confidence

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_confidence
```

````

````{py:method} get_proc_random_count()
:canonical: eon.akmcstate.AKMCState.get_proc_random_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_proc_random_count
```

````

````{py:method} inc_proc_random_count(procid)
:canonical: eon.akmcstate.AKMCState.inc_proc_random_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.inc_proc_random_count
```

````

````{py:method} reset_repeats()
:canonical: eon.akmcstate.AKMCState.reset_repeats

```{autodoc2-docstring} eon.akmcstate.AKMCState.reset_repeats
```

````

````{py:method} get_repeats()
:canonical: eon.akmcstate.AKMCState.get_repeats

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_repeats
```

````

````{py:method} inc_repeats()
:canonical: eon.akmcstate.AKMCState.inc_repeats

```{autodoc2-docstring} eon.akmcstate.AKMCState.inc_repeats
```

````

````{py:method} load_process_table()
:canonical: eon.akmcstate.AKMCState.load_process_table

```{autodoc2-docstring} eon.akmcstate.AKMCState.load_process_table
```

````

````{py:method} save_process_table()
:canonical: eon.akmcstate.AKMCState.save_process_table

```{autodoc2-docstring} eon.akmcstate.AKMCState.save_process_table
```

````

````{py:method} append_process_table(id, saddle_energy, prefactor, product, product_energy, product_prefactor, barrier, rate, repeats)
:canonical: eon.akmcstate.AKMCState.append_process_table

```{autodoc2-docstring} eon.akmcstate.AKMCState.append_process_table
```

````

````{py:method} update_lowest_barrier(barrier)
:canonical: eon.akmcstate.AKMCState.update_lowest_barrier

```{autodoc2-docstring} eon.akmcstate.AKMCState.update_lowest_barrier
```

````

````{py:method} get_lowest_barrier()
:canonical: eon.akmcstate.AKMCState.get_lowest_barrier

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_lowest_barrier
```

````

````{py:method} get_unique_saddle_count()
:canonical: eon.akmcstate.AKMCState.get_unique_saddle_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_unique_saddle_count
```

````

````{py:method} set_unique_saddle_count(num)
:canonical: eon.akmcstate.AKMCState.set_unique_saddle_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.set_unique_saddle_count
```

````

````{py:method} get_good_saddle_count()
:canonical: eon.akmcstate.AKMCState.get_good_saddle_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_good_saddle_count
```

````

````{py:method} set_good_saddle_count(num)
:canonical: eon.akmcstate.AKMCState.set_good_saddle_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.set_good_saddle_count
```

````

````{py:method} increment_time(dt, T_search)
:canonical: eon.akmcstate.AKMCState.increment_time

```{autodoc2-docstring} eon.akmcstate.AKMCState.increment_time
```

````

````{py:method} get_time()
:canonical: eon.akmcstate.AKMCState.get_time

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_time
```

````

````{py:method} get_time_by_temp()
:canonical: eon.akmcstate.AKMCState.get_time_by_temp

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_time_by_temp
```

````

````{py:method} get_number_of_searches()
:canonical: eon.akmcstate.AKMCState.get_number_of_searches

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_number_of_searches
```

````

````{py:method} get_total_saddle_count()
:canonical: eon.akmcstate.AKMCState.get_total_saddle_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_total_saddle_count
```

````

````{py:method} get_bad_saddle_count()
:canonical: eon.akmcstate.AKMCState.get_bad_saddle_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_bad_saddle_count
```

````

````{py:method} set_bad_saddle_count(num)
:canonical: eon.akmcstate.AKMCState.set_bad_saddle_count

```{autodoc2-docstring} eon.akmcstate.AKMCState.set_bad_saddle_count
```

````

````{py:method} register_bad_saddle(result, store=False, superbasin=None)
:canonical: eon.akmcstate.AKMCState.register_bad_saddle

```{autodoc2-docstring} eon.akmcstate.AKMCState.register_bad_saddle
```

````

````{py:method} get_process_reactant(id)
:canonical: eon.akmcstate.AKMCState.get_process_reactant

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_process_reactant
```

````

````{py:method} get_process_saddle(id)
:canonical: eon.akmcstate.AKMCState.get_process_saddle

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_process_saddle
```

````

````{py:method} get_process_product(id)
:canonical: eon.akmcstate.AKMCState.get_process_product

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_process_product
```

````

````{py:method} get_process_mode(id)
:canonical: eon.akmcstate.AKMCState.get_process_mode

```{autodoc2-docstring} eon.akmcstate.AKMCState.get_process_mode
```

````

````{py:method} proc_saddle_path(id)
:canonical: eon.akmcstate.AKMCState.proc_saddle_path

```{autodoc2-docstring} eon.akmcstate.AKMCState.proc_saddle_path
```

````

````{py:method} proc_mode_path(id)
:canonical: eon.akmcstate.AKMCState.proc_mode_path

```{autodoc2-docstring} eon.akmcstate.AKMCState.proc_mode_path
```

````

`````

````{py:function} lambertw(z)
:canonical: eon.akmcstate.lambertw

```{autodoc2-docstring} eon.akmcstate.lambertw
```
````
