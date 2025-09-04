# {py:mod}`eon.askmc`

```{py:module} eon.askmc
```

```{autodoc2-docstring} eon.askmc
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`ASKMC <eon.askmc.ASKMC>`
  - ```{autodoc2-docstring} eon.askmc.ASKMC
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.askmc.logger>`
  - ```{autodoc2-docstring} eon.askmc.logger
    :summary:
    ```
* - {py:obj}`processtable_head_fmt <eon.askmc.processtable_head_fmt>`
  - ```{autodoc2-docstring} eon.askmc.processtable_head_fmt
    :summary:
    ```
* - {py:obj}`mod_processtable_head_fmt <eon.askmc.mod_processtable_head_fmt>`
  - ```{autodoc2-docstring} eon.askmc.mod_processtable_head_fmt
    :summary:
    ```
* - {py:obj}`processtable_header <eon.askmc.processtable_header>`
  - ```{autodoc2-docstring} eon.askmc.processtable_header
    :summary:
    ```
* - {py:obj}`mod_processtable_header <eon.askmc.mod_processtable_header>`
  - ```{autodoc2-docstring} eon.askmc.mod_processtable_header
    :summary:
    ```
* - {py:obj}`processtable_line <eon.askmc.processtable_line>`
  - ```{autodoc2-docstring} eon.askmc.processtable_line
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.askmc.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.askmc.logger
```

````

````{py:data} processtable_head_fmt
:canonical: eon.askmc.processtable_head_fmt
:value: >
   '%7s %16s %11s %9s %16s %17s %8s %12s %7s\n'

```{autodoc2-docstring} eon.askmc.processtable_head_fmt
```

````

````{py:data} mod_processtable_head_fmt
:canonical: eon.askmc.mod_processtable_head_fmt
:value: >
   '%7s %16s %11s %9s %16s %17s %8s %12s %12s\n'

```{autodoc2-docstring} eon.askmc.mod_processtable_head_fmt
```

````

````{py:data} processtable_header
:canonical: eon.askmc.processtable_header
:value: >
   None

```{autodoc2-docstring} eon.askmc.processtable_header
```

````

````{py:data} mod_processtable_header
:canonical: eon.askmc.mod_processtable_header
:value: >
   None

```{autodoc2-docstring} eon.askmc.mod_processtable_header
```

````

````{py:data} processtable_line
:canonical: eon.askmc.processtable_line
:value: >
   '%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n'

```{autodoc2-docstring} eon.askmc.processtable_line
```

````

`````{py:class} ASKMC(kT, states, confidence, alpha, gamma, barrier_test_on, connection_test_on, sb_recycling_on, path_root, thermal_window, recycle_path)
:canonical: eon.askmc.ASKMC

```{autodoc2-docstring} eon.askmc.ASKMC
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.askmc.ASKMC.__init__
```

````{py:method} compile_process_table(current_state)
:canonical: eon.askmc.ASKMC.compile_process_table

```{autodoc2-docstring} eon.askmc.ASKMC.compile_process_table
```

````

````{py:method} get_ratetable(current_state)
:canonical: eon.askmc.ASKMC.get_ratetable

```{autodoc2-docstring} eon.askmc.ASKMC.get_ratetable
```

````

````{py:method} get_askmc_metadata()
:canonical: eon.askmc.ASKMC.get_askmc_metadata

```{autodoc2-docstring} eon.askmc.ASKMC.get_askmc_metadata
```

````

````{py:method} save_askmc_metadata(sb_check_count, num_rate_changes)
:canonical: eon.askmc.ASKMC.save_askmc_metadata

```{autodoc2-docstring} eon.askmc.ASKMC.save_askmc_metadata
```

````

````{py:method} get_real_process_table(current_state)
:canonical: eon.askmc.ASKMC.get_real_process_table

```{autodoc2-docstring} eon.askmc.ASKMC.get_real_process_table
```

````

````{py:method} get_modified_process_table(current_state)
:canonical: eon.askmc.ASKMC.get_modified_process_table

```{autodoc2-docstring} eon.askmc.ASKMC.get_modified_process_table
```

````

````{py:method} save_modified_process_table(current_state, current_state_mod_procs)
:canonical: eon.askmc.ASKMC.save_modified_process_table

```{autodoc2-docstring} eon.askmc.ASKMC.save_modified_process_table
```

````

````{py:method} append_modified_process_table(current_state, process_id, saddle_energy, prefactor, product, product_energy, product_prefactor, barrier, rate, view_count)
:canonical: eon.askmc.ASKMC.append_modified_process_table

```{autodoc2-docstring} eon.askmc.ASKMC.append_modified_process_table
```

````

````{py:method} get_process_id(current_state_procs, next_state_num, flag)
:canonical: eon.askmc.ASKMC.get_process_id

```{autodoc2-docstring} eon.askmc.ASKMC.get_process_id
```

````

````{py:method} register_transition(current_state, next_state)
:canonical: eon.askmc.ASKMC.register_transition

```{autodoc2-docstring} eon.askmc.ASKMC.register_transition
```

````

````{py:method} in_array(test, array)
:canonical: eon.askmc.ASKMC.in_array

```{autodoc2-docstring} eon.askmc.ASKMC.in_array
```

````

````{py:method} is_equal(a, b)
:canonical: eon.askmc.ASKMC.is_equal

```{autodoc2-docstring} eon.askmc.ASKMC.is_equal
```

````

````{py:method} edgelist_to_statelist()
:canonical: eon.askmc.ASKMC.edgelist_to_statelist

```{autodoc2-docstring} eon.askmc.ASKMC.edgelist_to_statelist
```

````

````{py:method} missed_connections(origEtrans)
:canonical: eon.askmc.ASKMC.missed_connections

```{autodoc2-docstring} eon.askmc.ASKMC.missed_connections
```

````

````{py:method} missed_low_barriers(origEtrans)
:canonical: eon.askmc.ASKMC.missed_low_barriers

```{autodoc2-docstring} eon.askmc.ASKMC.missed_low_barriers
```

````

````{py:method} locsearch(current_state, origEtrans)
:canonical: eon.askmc.ASKMC.locsearch

```{autodoc2-docstring} eon.askmc.ASKMC.locsearch
```

````

````{py:method} raiseup(current_state, next_state, sb_check_count, num_rate_changes)
:canonical: eon.askmc.ASKMC.raiseup

```{autodoc2-docstring} eon.askmc.ASKMC.raiseup
```

````

`````
