# {py:mod}`eon.akmcstatelist`

```{py:module} eon.akmcstatelist
```

```{autodoc2-docstring} eon.akmcstatelist
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`AKMCStateList <eon.akmcstatelist.AKMCStateList>`
  - ```{autodoc2-docstring} eon.akmcstatelist.AKMCStateList
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.akmcstatelist.logger>`
  - ```{autodoc2-docstring} eon.akmcstatelist.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.akmcstatelist.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.akmcstatelist.logger
```

````

`````{py:class} AKMCStateList(kT, thermal_window, max_thermal_window, initial_state=None, filter_hole=False, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.akmcstatelist.AKMCStateList

Bases: {py:obj}`eon.statelist.StateList`

```{autodoc2-docstring} eon.akmcstatelist.AKMCStateList
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.akmcstatelist.AKMCStateList.__init__
```

````{py:method} register_process(reactant_number, product_number, process_id)
:canonical: eon.akmcstatelist.AKMCStateList.register_process

```{autodoc2-docstring} eon.akmcstatelist.AKMCStateList.register_process
```

````

````{py:method} connect_states(states)
:canonical: eon.akmcstatelist.AKMCStateList.connect_states

```{autodoc2-docstring} eon.akmcstatelist.AKMCStateList.connect_states
```

````

````{py:method} connect_state_sets(states1, states2)
:canonical: eon.akmcstatelist.AKMCStateList.connect_state_sets

```{autodoc2-docstring} eon.akmcstatelist.AKMCStateList.connect_state_sets
```

````

`````
