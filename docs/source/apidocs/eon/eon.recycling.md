# {py:mod}`eon.recycling`

```{py:module} eon.recycling
```

```{autodoc2-docstring} eon.recycling
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`SB_Recycling <eon.recycling.SB_Recycling>`
  - ```{autodoc2-docstring} eon.recycling.SB_Recycling
    :summary:
    ```
* - {py:obj}`Recycling <eon.recycling.Recycling>`
  - ```{autodoc2-docstring} eon.recycling.Recycling
    :summary:
    ```
````

### API

`````{py:class} SB_Recycling(states, previous_state, current_state, move_distance, recycle_save, path, sb_scheme, superbasining)
:canonical: eon.recycling.SB_Recycling

```{autodoc2-docstring} eon.recycling.SB_Recycling
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.recycling.SB_Recycling.__init__
```

````{py:method} get_process_id(current_state_procs, next_state_num)
:canonical: eon.recycling.SB_Recycling.get_process_id

```{autodoc2-docstring} eon.recycling.SB_Recycling.get_process_id
```

````

````{py:method} make_suggestion()
:canonical: eon.recycling.SB_Recycling.make_suggestion

```{autodoc2-docstring} eon.recycling.SB_Recycling.make_suggestion
```

````

````{py:method} write_metadata()
:canonical: eon.recycling.SB_Recycling.write_metadata

```{autodoc2-docstring} eon.recycling.SB_Recycling.write_metadata
```

````

````{py:method} read_metadata()
:canonical: eon.recycling.SB_Recycling.read_metadata

```{autodoc2-docstring} eon.recycling.SB_Recycling.read_metadata
```

````

````{py:method} generate_corresponding_states()
:canonical: eon.recycling.SB_Recycling.generate_corresponding_states

```{autodoc2-docstring} eon.recycling.SB_Recycling.generate_corresponding_states
```

````

`````

`````{py:class} Recycling(states, suggested_ref_state, new_state, move_distance, save=False, from_sb=False, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.recycling.Recycling

```{autodoc2-docstring} eon.recycling.Recycling
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.recycling.Recycling.__init__
```

````{py:method} make_suggestion()
:canonical: eon.recycling.Recycling.make_suggestion

```{autodoc2-docstring} eon.recycling.Recycling.make_suggestion
```

````

````{py:method} read_recycling_metadata()
:canonical: eon.recycling.Recycling.read_recycling_metadata

```{autodoc2-docstring} eon.recycling.Recycling.read_recycling_metadata
```

````

````{py:method} write_recycling_metadata()
:canonical: eon.recycling.Recycling.write_recycling_metadata

```{autodoc2-docstring} eon.recycling.Recycling.write_recycling_metadata
```

````

````{py:method} get_moved_indices()
:canonical: eon.recycling.Recycling.get_moved_indices

```{autodoc2-docstring} eon.recycling.Recycling.get_moved_indices
```

````

`````
