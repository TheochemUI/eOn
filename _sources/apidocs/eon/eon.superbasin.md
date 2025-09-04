# {py:mod}`eon.superbasin`

```{py:module} eon.superbasin
```

```{autodoc2-docstring} eon.superbasin
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Superbasin <eon.superbasin.Superbasin>`
  - ```{autodoc2-docstring} eon.superbasin.Superbasin
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.superbasin.logger>`
  - ```{autodoc2-docstring} eon.superbasin.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.superbasin.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.superbasin.logger
```

````

`````{py:class} Superbasin(path, id, state_list=None, get_state=None)
:canonical: eon.superbasin.Superbasin

```{autodoc2-docstring} eon.superbasin.Superbasin
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.superbasin.Superbasin.__init__
```

````{py:method} step(entry_state, get_product_state)
:canonical: eon.superbasin.Superbasin.step

```{autodoc2-docstring} eon.superbasin.Superbasin.step
```

````

````{py:method} contains_state(state)
:canonical: eon.superbasin.Superbasin.contains_state

```{autodoc2-docstring} eon.superbasin.Superbasin.contains_state
```

````

````{py:method} write_data()
:canonical: eon.superbasin.Superbasin.write_data

```{autodoc2-docstring} eon.superbasin.Superbasin.write_data
```

````

````{py:method} read_data(get_state)
:canonical: eon.superbasin.Superbasin.read_data

```{autodoc2-docstring} eon.superbasin.Superbasin.read_data
```

````

````{py:method} delete(storage=None)
:canonical: eon.superbasin.Superbasin.delete

```{autodoc2-docstring} eon.superbasin.Superbasin.delete
```

````

````{py:method} get_confidence()
:canonical: eon.superbasin.Superbasin.get_confidence

```{autodoc2-docstring} eon.superbasin.Superbasin.get_confidence
```

````

````{py:method} get_lowest_confidence_state()
:canonical: eon.superbasin.Superbasin.get_lowest_confidence_state

```{autodoc2-docstring} eon.superbasin.Superbasin.get_lowest_confidence_state
```

````

`````
