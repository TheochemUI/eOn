# {py:mod}`eon.explorer`

```{py:module} eon.explorer
```

```{autodoc2-docstring} eon.explorer
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Explorer <eon.explorer.Explorer>`
  - ```{autodoc2-docstring} eon.explorer.Explorer
    :summary:
    ```
* - {py:obj}`MinModeExplorer <eon.explorer.MinModeExplorer>`
  - ```{autodoc2-docstring} eon.explorer.MinModeExplorer
    :summary:
    ```
* - {py:obj}`ClientMinModeExplorer <eon.explorer.ClientMinModeExplorer>`
  - ```{autodoc2-docstring} eon.explorer.ClientMinModeExplorer
    :summary:
    ```
* - {py:obj}`ServerMinModeExplorer <eon.explorer.ServerMinModeExplorer>`
  - ```{autodoc2-docstring} eon.explorer.ServerMinModeExplorer
    :summary:
    ```
* - {py:obj}`ProcessSearch <eon.explorer.ProcessSearch>`
  - ```{autodoc2-docstring} eon.explorer.ProcessSearch
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`get_minmodexplorer <eon.explorer.get_minmodexplorer>`
  - ```{autodoc2-docstring} eon.explorer.get_minmodexplorer
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.explorer.logger>`
  - ```{autodoc2-docstring} eon.explorer.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.explorer.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.explorer.logger
```

````

````{py:function} get_minmodexplorer(config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.explorer.get_minmodexplorer

```{autodoc2-docstring} eon.explorer.get_minmodexplorer
```
````

`````{py:class} Explorer(superbasin=None, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.explorer.Explorer

```{autodoc2-docstring} eon.explorer.Explorer
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.explorer.Explorer.__init__
```

````{py:method} load_wuid()
:canonical: eon.explorer.Explorer.load_wuid

```{autodoc2-docstring} eon.explorer.Explorer.load_wuid
```

````

````{py:method} save_wuid()
:canonical: eon.explorer.Explorer.save_wuid

```{autodoc2-docstring} eon.explorer.Explorer.save_wuid
```

````

`````

`````{py:class} MinModeExplorer(states, previous_state, state, superbasin=None, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.explorer.MinModeExplorer

Bases: {py:obj}`eon.explorer.Explorer`

```{autodoc2-docstring} eon.explorer.MinModeExplorer
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.explorer.MinModeExplorer.__init__
```

````{py:method} explore()
:canonical: eon.explorer.MinModeExplorer.explore

```{autodoc2-docstring} eon.explorer.MinModeExplorer.explore
```

````

````{py:method} generate_displacement()
:canonical: eon.explorer.MinModeExplorer.generate_displacement

```{autodoc2-docstring} eon.explorer.MinModeExplorer.generate_displacement
```

````

`````

`````{py:class} ClientMinModeExplorer(states, previous_state, state, superbasin=None)
:canonical: eon.explorer.ClientMinModeExplorer

Bases: {py:obj}`eon.explorer.MinModeExplorer`

```{autodoc2-docstring} eon.explorer.ClientMinModeExplorer
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.explorer.ClientMinModeExplorer.__init__
```

````{py:method} make_jobs()
:canonical: eon.explorer.ClientMinModeExplorer.make_jobs

```{autodoc2-docstring} eon.explorer.ClientMinModeExplorer.make_jobs
```

````

````{py:method} register_results()
:canonical: eon.explorer.ClientMinModeExplorer.register_results

```{autodoc2-docstring} eon.explorer.ClientMinModeExplorer.register_results
```

````

`````

`````{py:class} ServerMinModeExplorer(states, previous_state, state, superbasin=None)
:canonical: eon.explorer.ServerMinModeExplorer

Bases: {py:obj}`eon.explorer.MinModeExplorer`

```{autodoc2-docstring} eon.explorer.ServerMinModeExplorer
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.explorer.ServerMinModeExplorer.__init__
```

````{py:method} save()
:canonical: eon.explorer.ServerMinModeExplorer.save

```{autodoc2-docstring} eon.explorer.ServerMinModeExplorer.save
```

````

````{py:method} explore()
:canonical: eon.explorer.ServerMinModeExplorer.explore

```{autodoc2-docstring} eon.explorer.ServerMinModeExplorer.explore
```

````

````{py:method} register_results()
:canonical: eon.explorer.ServerMinModeExplorer.register_results

```{autodoc2-docstring} eon.explorer.ServerMinModeExplorer.register_results
```

````

````{py:method} make_jobs()
:canonical: eon.explorer.ServerMinModeExplorer.make_jobs

```{autodoc2-docstring} eon.explorer.ServerMinModeExplorer.make_jobs
```

````

`````

`````{py:class} ProcessSearch(reactant, displacement, mode, disp_type, search_id, state_number)
:canonical: eon.explorer.ProcessSearch

```{autodoc2-docstring} eon.explorer.ProcessSearch
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.explorer.ProcessSearch.__init__
```

````{py:method} get_job(state_number)
:canonical: eon.explorer.ProcessSearch.get_job

```{autodoc2-docstring} eon.explorer.ProcessSearch.get_job
```

````

````{py:method} process_result(result)
:canonical: eon.explorer.ProcessSearch.process_result

```{autodoc2-docstring} eon.explorer.ProcessSearch.process_result
```

````

````{py:method} build_result()
:canonical: eon.explorer.ProcessSearch.build_result

```{autodoc2-docstring} eon.explorer.ProcessSearch.build_result
```

````

````{py:method} get_saddle()
:canonical: eon.explorer.ProcessSearch.get_saddle

```{autodoc2-docstring} eon.explorer.ProcessSearch.get_saddle
```

````

````{py:method} get_saddle_file()
:canonical: eon.explorer.ProcessSearch.get_saddle_file

```{autodoc2-docstring} eon.explorer.ProcessSearch.get_saddle_file
```

````

````{py:method} start_minimization(which_min)
:canonical: eon.explorer.ProcessSearch.start_minimization

```{autodoc2-docstring} eon.explorer.ProcessSearch.start_minimization
```

````

````{py:method} finish_minimization(result)
:canonical: eon.explorer.ProcessSearch.finish_minimization

```{autodoc2-docstring} eon.explorer.ProcessSearch.finish_minimization
```

````

````{py:method} start_search()
:canonical: eon.explorer.ProcessSearch.start_search

```{autodoc2-docstring} eon.explorer.ProcessSearch.start_search
```

````

````{py:method} finish_search(result)
:canonical: eon.explorer.ProcessSearch.finish_search

```{autodoc2-docstring} eon.explorer.ProcessSearch.finish_search
```

````

````{py:method} save_result(result)
:canonical: eon.explorer.ProcessSearch.save_result

```{autodoc2-docstring} eon.explorer.ProcessSearch.save_result
```

````

````{py:method} load_result(result_name)
:canonical: eon.explorer.ProcessSearch.load_result

```{autodoc2-docstring} eon.explorer.ProcessSearch.load_result
```

````

`````
