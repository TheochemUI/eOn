# {py:mod}`eon.superbasinscheme`

```{py:module} eon.superbasinscheme
```

```{autodoc2-docstring} eon.superbasinscheme
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`SuperbasinScheme <eon.superbasinscheme.SuperbasinScheme>`
  - ```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme
    :summary:
    ```
* - {py:obj}`TransitionCounting <eon.superbasinscheme.TransitionCounting>`
  - ```{autodoc2-docstring} eon.superbasinscheme.TransitionCounting
    :summary:
    ```
* - {py:obj}`EnergyLevel <eon.superbasinscheme.EnergyLevel>`
  -
* - {py:obj}`RateThreshold <eon.superbasinscheme.RateThreshold>`
  - ```{autodoc2-docstring} eon.superbasinscheme.RateThreshold
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.superbasinscheme.logger>`
  - ```{autodoc2-docstring} eon.superbasinscheme.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.superbasinscheme.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.superbasinscheme.logger
```

````

`````{py:class} SuperbasinScheme(superbasin_path, states, kT)
:canonical: eon.superbasinscheme.SuperbasinScheme

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme.__init__
```

````{py:method} get_containing_superbasin(state)
:canonical: eon.superbasinscheme.SuperbasinScheme.get_containing_superbasin

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme.get_containing_superbasin
```

````

````{py:method} make_basin_from_sets(start_state, end_state)
:canonical: eon.superbasinscheme.SuperbasinScheme.make_basin_from_sets

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme.make_basin_from_sets
```

````

````{py:method} make_basin(merge_states)
:canonical: eon.superbasinscheme.SuperbasinScheme.make_basin

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme.make_basin
```

````

````{py:method} register_transition(start_state, end_state)
:canonical: eon.superbasinscheme.SuperbasinScheme.register_transition
:abstractmethod:

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme.register_transition
```

````

````{py:method} write_data()
:canonical: eon.superbasinscheme.SuperbasinScheme.write_data
:abstractmethod:

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme.write_data
```

````

````{py:method} read_data()
:canonical: eon.superbasinscheme.SuperbasinScheme.read_data
:abstractmethod:

```{autodoc2-docstring} eon.superbasinscheme.SuperbasinScheme.read_data
```

````

`````

`````{py:class} TransitionCounting(superbasin_path, states, kT, num_transitions)
:canonical: eon.superbasinscheme.TransitionCounting

Bases: {py:obj}`eon.superbasinscheme.SuperbasinScheme`

```{autodoc2-docstring} eon.superbasinscheme.TransitionCounting
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.superbasinscheme.TransitionCounting.__init__
```

````{py:method} register_transition(start_state, end_state)
:canonical: eon.superbasinscheme.TransitionCounting.register_transition

```{autodoc2-docstring} eon.superbasinscheme.TransitionCounting.register_transition
```

````

````{py:method} write_data()
:canonical: eon.superbasinscheme.TransitionCounting.write_data

```{autodoc2-docstring} eon.superbasinscheme.TransitionCounting.write_data
```

````

````{py:method} read_data()
:canonical: eon.superbasinscheme.TransitionCounting.read_data

```{autodoc2-docstring} eon.superbasinscheme.TransitionCounting.read_data
```

````

````{py:method} get_count(state)
:canonical: eon.superbasinscheme.TransitionCounting.get_count

```{autodoc2-docstring} eon.superbasinscheme.TransitionCounting.get_count
```

````

`````

`````{py:class} EnergyLevel(superbasin_path, states, kT, energy_increment)
:canonical: eon.superbasinscheme.EnergyLevel

Bases: {py:obj}`eon.superbasinscheme.SuperbasinScheme`

````{py:method} get_energy_increment(e_min, e_global_min, e_saddle)
:canonical: eon.superbasinscheme.EnergyLevel.get_energy_increment

```{autodoc2-docstring} eon.superbasinscheme.EnergyLevel.get_energy_increment
```

````

````{py:method} get_statelist(state)
:canonical: eon.superbasinscheme.EnergyLevel.get_statelist

```{autodoc2-docstring} eon.superbasinscheme.EnergyLevel.get_statelist
```

````

````{py:method} register_transition(start_state, end_state)
:canonical: eon.superbasinscheme.EnergyLevel.register_transition

```{autodoc2-docstring} eon.superbasinscheme.EnergyLevel.register_transition
```

````

````{py:method} read_data()
:canonical: eon.superbasinscheme.EnergyLevel.read_data

```{autodoc2-docstring} eon.superbasinscheme.EnergyLevel.read_data
```

````

````{py:method} write_data()
:canonical: eon.superbasinscheme.EnergyLevel.write_data

```{autodoc2-docstring} eon.superbasinscheme.EnergyLevel.write_data
```

````

`````

`````{py:class} RateThreshold(superbasin_path, states, kT, rate_threshold)
:canonical: eon.superbasinscheme.RateThreshold

Bases: {py:obj}`eon.superbasinscheme.SuperbasinScheme`

```{autodoc2-docstring} eon.superbasinscheme.RateThreshold
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.superbasinscheme.RateThreshold.__init__
```

````{py:method} register_transition(start_state, end_state)
:canonical: eon.superbasinscheme.RateThreshold.register_transition

```{autodoc2-docstring} eon.superbasinscheme.RateThreshold.register_transition
```

````

````{py:method} write_data()
:canonical: eon.superbasinscheme.RateThreshold.write_data

```{autodoc2-docstring} eon.superbasinscheme.RateThreshold.write_data
```

````

````{py:method} read_data()
:canonical: eon.superbasinscheme.RateThreshold.read_data

```{autodoc2-docstring} eon.superbasinscheme.RateThreshold.read_data
```

````

`````
