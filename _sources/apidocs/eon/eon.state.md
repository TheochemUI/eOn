# {py:mod}`eon.state`

```{py:module} eon.state
```

```{autodoc2-docstring} eon.state
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`State <eon.state.State>`
  - ```{autodoc2-docstring} eon.state.State
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.state.logger>`
  - ```{autodoc2-docstring} eon.state.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.state.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.state.logger
```

````

`````{py:class} State(statepath, statenumber, statelist, previous_state_num=-1, reactant_path=None, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.state.State

```{autodoc2-docstring} eon.state.State
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.state.State.__init__
```

````{py:method} add_process(result)
:canonical: eon.state.State.add_process

```{autodoc2-docstring} eon.state.State.add_process
```

````

````{py:method} get_energy()
:canonical: eon.state.State.get_energy

```{autodoc2-docstring} eon.state.State.get_energy
```

````

````{py:method} set_energy(e)
:canonical: eon.state.State.set_energy

```{autodoc2-docstring} eon.state.State.set_energy
```

````

````{py:method} get_reactant()
:canonical: eon.state.State.get_reactant

```{autodoc2-docstring} eon.state.State.get_reactant
```

````

````{py:method} get_process_reactant(id)
:canonical: eon.state.State.get_process_reactant

```{autodoc2-docstring} eon.state.State.get_process_reactant
```

````

````{py:method} get_process_product(id)
:canonical: eon.state.State.get_process_product

```{autodoc2-docstring} eon.state.State.get_process_product
```

````

````{py:method} proc_reactant_path(id)
:canonical: eon.state.State.proc_reactant_path

```{autodoc2-docstring} eon.state.State.proc_reactant_path
```

````

````{py:method} proc_product_path(id)
:canonical: eon.state.State.proc_product_path

```{autodoc2-docstring} eon.state.State.proc_product_path
```

````

````{py:method} proc_results_path(id)
:canonical: eon.state.State.proc_results_path

```{autodoc2-docstring} eon.state.State.proc_results_path
```

````

````{py:method} proc_stdout_path(id)
:canonical: eon.state.State.proc_stdout_path

```{autodoc2-docstring} eon.state.State.proc_stdout_path
```

````

````{py:method} tar_procdata()
:canonical: eon.state.State.tar_procdata

```{autodoc2-docstring} eon.state.State.tar_procdata
```

````

````{py:method} get_process(id)
:canonical: eon.state.State.get_process

```{autodoc2-docstring} eon.state.State.get_process
```

````

````{py:method} get_process_ids()
:canonical: eon.state.State.get_process_ids

```{autodoc2-docstring} eon.state.State.get_process_ids
```

````

````{py:method} get_previous_state()
:canonical: eon.state.State.get_previous_state

```{autodoc2-docstring} eon.state.State.get_previous_state
```

````

````{py:method} get_num_procs()
:canonical: eon.state.State.get_num_procs

```{autodoc2-docstring} eon.state.State.get_num_procs
```

````

````{py:method} get_process_table()
:canonical: eon.state.State.get_process_table

```{autodoc2-docstring} eon.state.State.get_process_table
```

````

`````
