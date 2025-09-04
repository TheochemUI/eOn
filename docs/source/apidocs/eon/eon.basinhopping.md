# {py:mod}`eon.basinhopping`

```{py:module} eon.basinhopping
```

```{autodoc2-docstring} eon.basinhopping
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`BHStates <eon.basinhopping.BHStates>`
  - ```{autodoc2-docstring} eon.basinhopping.BHStates
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`basinhopping <eon.basinhopping.basinhopping>`
  - ```{autodoc2-docstring} eon.basinhopping.basinhopping
    :summary:
    ```
* - {py:obj}`make_searches <eon.basinhopping.make_searches>`
  - ```{autodoc2-docstring} eon.basinhopping.make_searches
    :summary:
    ```
* - {py:obj}`register_results <eon.basinhopping.register_results>`
  - ```{autodoc2-docstring} eon.basinhopping.register_results
    :summary:
    ```
* - {py:obj}`main <eon.basinhopping.main>`
  - ```{autodoc2-docstring} eon.basinhopping.main
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.basinhopping.logger>`
  - ```{autodoc2-docstring} eon.basinhopping.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.basinhopping.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.basinhopping.logger
```

````

`````{py:class} BHStates()
:canonical: eon.basinhopping.BHStates

```{autodoc2-docstring} eon.basinhopping.BHStates
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.basinhopping.BHStates.__init__
```

````{py:method} get_random_minimum()
:canonical: eon.basinhopping.BHStates.get_random_minimum

```{autodoc2-docstring} eon.basinhopping.BHStates.get_random_minimum
```

````

````{py:method} add_state(result_files, result_info)
:canonical: eon.basinhopping.BHStates.add_state

```{autodoc2-docstring} eon.basinhopping.BHStates.add_state
```

````

`````

````{py:function} basinhopping()
:canonical: eon.basinhopping.basinhopping

```{autodoc2-docstring} eon.basinhopping.basinhopping
```
````

````{py:function} make_searches(comm, wuid, bhstates)
:canonical: eon.basinhopping.make_searches

```{autodoc2-docstring} eon.basinhopping.make_searches
```
````

````{py:function} register_results(comm, bhstates)
:canonical: eon.basinhopping.register_results

```{autodoc2-docstring} eon.basinhopping.register_results
```
````

````{py:function} main()
:canonical: eon.basinhopping.main

```{autodoc2-docstring} eon.basinhopping.main
```
````
