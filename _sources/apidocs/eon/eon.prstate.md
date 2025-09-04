# {py:mod}`eon.prstate`

```{py:module} eon.prstate
```

```{autodoc2-docstring} eon.prstate
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`PRState <eon.prstate.PRState>`
  - ```{autodoc2-docstring} eon.prstate.PRState
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.prstate.logger>`
  - ```{autodoc2-docstring} eon.prstate.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.prstate.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.prstate.logger
```

````

`````{py:class} PRState(statepath, statenumber, statelist, previous_state_num=-1, reactant_path=None)
:canonical: eon.prstate.PRState

Bases: {py:obj}`eon.state.State`

```{autodoc2-docstring} eon.prstate.PRState
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.prstate.PRState.__init__
```

````{py:attribute} processtable_head_fmt
:canonical: eon.prstate.PRState.processtable_head_fmt
:value: >
   '%7s %9s %16s %12s\n'

```{autodoc2-docstring} eon.prstate.PRState.processtable_head_fmt
```

````

````{py:attribute} processtable_header
:canonical: eon.prstate.PRState.processtable_header
:value: >
   None

```{autodoc2-docstring} eon.prstate.PRState.processtable_header
```

````

````{py:attribute} processtable_line
:canonical: eon.prstate.PRState.processtable_line
:value: >
   '%7d %9d %16.5f %12.5e\n'

```{autodoc2-docstring} eon.prstate.PRState.processtable_line
```

````

````{py:attribute} search_result_header
:canonical: eon.prstate.PRState.search_result_header
:value: >
   None

```{autodoc2-docstring} eon.prstate.PRState.search_result_header
```

````

````{py:method} add_process(result)
:canonical: eon.prstate.PRState.add_process

```{autodoc2-docstring} eon.prstate.PRState.add_process
```

````

````{py:method} load_process_table()
:canonical: eon.prstate.PRState.load_process_table

```{autodoc2-docstring} eon.prstate.PRState.load_process_table
```

````

````{py:method} save_process_table()
:canonical: eon.prstate.PRState.save_process_table

```{autodoc2-docstring} eon.prstate.PRState.save_process_table
```

````

````{py:method} append_process_table(id, product, product_energy, time)
:canonical: eon.prstate.PRState.append_process_table

```{autodoc2-docstring} eon.prstate.PRState.append_process_table
```

````

````{py:method} get_time()
:canonical: eon.prstate.PRState.get_time

```{autodoc2-docstring} eon.prstate.PRState.get_time
```

````

````{py:method} inc_time(timeinc)
:canonical: eon.prstate.PRState.inc_time

```{autodoc2-docstring} eon.prstate.PRState.inc_time
```

````

````{py:method} zero_time()
:canonical: eon.prstate.PRState.zero_time

```{autodoc2-docstring} eon.prstate.PRState.zero_time
```

````

`````
