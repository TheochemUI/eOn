# {py:mod}`eon.statelist`

```{py:module} eon.statelist
```

```{autodoc2-docstring} eon.statelist
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`StateList <eon.statelist.StateList>`
  - ```{autodoc2-docstring} eon.statelist.StateList
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.statelist.logger>`
  - ```{autodoc2-docstring} eon.statelist.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.statelist.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.statelist.logger
```

````

`````{py:class} StateList(StateClass, initial_state=None, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.statelist.StateList

```{autodoc2-docstring} eon.statelist.StateList
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.statelist.StateList.__init__
```

````{py:method} get_num_states()
:canonical: eon.statelist.StateList.get_num_states

```{autodoc2-docstring} eon.statelist.StateList.get_num_states
```

````

````{py:method} get_product_state(state_number, process_id)
:canonical: eon.statelist.StateList.get_product_state

```{autodoc2-docstring} eon.statelist.StateList.get_product_state
```

````

````{py:method} get_state(state_number)
:canonical: eon.statelist.StateList.get_state

```{autodoc2-docstring} eon.statelist.StateList.get_state
```

````

````{py:method} load_state_table()
:canonical: eon.statelist.StateList.load_state_table

```{autodoc2-docstring} eon.statelist.StateList.load_state_table
```

````

````{py:method} save_state_table()
:canonical: eon.statelist.StateList.save_state_table

```{autodoc2-docstring} eon.statelist.StateList.save_state_table
```

````

````{py:method} append_state_table(energy)
:canonical: eon.statelist.StateList.append_state_table

```{autodoc2-docstring} eon.statelist.StateList.append_state_table
```

````

````{py:method} state_path(state_number)
:canonical: eon.statelist.StateList.state_path

```{autodoc2-docstring} eon.statelist.StateList.state_path
```

````

`````
