# {py:mod}`eon.escaperate`

```{py:module} eon.escaperate
```

```{autodoc2-docstring} eon.escaperate
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`parallelreplica <eon.escaperate.parallelreplica>`
  - ```{autodoc2-docstring} eon.escaperate.parallelreplica
    :summary:
    ```
* - {py:obj}`step <eon.escaperate.step>`
  - ```{autodoc2-docstring} eon.escaperate.step
    :summary:
    ```
* - {py:obj}`get_statelist <eon.escaperate.get_statelist>`
  - ```{autodoc2-docstring} eon.escaperate.get_statelist
    :summary:
    ```
* - {py:obj}`get_pr_metadata <eon.escaperate.get_pr_metadata>`
  - ```{autodoc2-docstring} eon.escaperate.get_pr_metadata
    :summary:
    ```
* - {py:obj}`write_pr_metadata <eon.escaperate.write_pr_metadata>`
  - ```{autodoc2-docstring} eon.escaperate.write_pr_metadata
    :summary:
    ```
* - {py:obj}`make_searches <eon.escaperate.make_searches>`
  - ```{autodoc2-docstring} eon.escaperate.make_searches
    :summary:
    ```
* - {py:obj}`register_results <eon.escaperate.register_results>`
  - ```{autodoc2-docstring} eon.escaperate.register_results
    :summary:
    ```
* - {py:obj}`main <eon.escaperate.main>`
  - ```{autodoc2-docstring} eon.escaperate.main
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.escaperate.logger>`
  - ```{autodoc2-docstring} eon.escaperate.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.escaperate.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.escaperate.logger
```

````

````{py:function} parallelreplica()
:canonical: eon.escaperate.parallelreplica

```{autodoc2-docstring} eon.escaperate.parallelreplica
```
````

````{py:function} step(current_time, current_state, states, transition)
:canonical: eon.escaperate.step

```{autodoc2-docstring} eon.escaperate.step
```
````

````{py:function} get_statelist()
:canonical: eon.escaperate.get_statelist

```{autodoc2-docstring} eon.escaperate.get_statelist
```
````

````{py:function} get_pr_metadata()
:canonical: eon.escaperate.get_pr_metadata

```{autodoc2-docstring} eon.escaperate.get_pr_metadata
```
````

````{py:function} write_pr_metadata(parser, current_state_num, time, wuid)
:canonical: eon.escaperate.write_pr_metadata

```{autodoc2-docstring} eon.escaperate.write_pr_metadata
```
````

````{py:function} make_searches(comm, current_state, wuid)
:canonical: eon.escaperate.make_searches

```{autodoc2-docstring} eon.escaperate.make_searches
```
````

````{py:function} register_results(comm, current_state, states)
:canonical: eon.escaperate.register_results

```{autodoc2-docstring} eon.escaperate.register_results
```
````

````{py:function} main()
:canonical: eon.escaperate.main

```{autodoc2-docstring} eon.escaperate.main
```
````
