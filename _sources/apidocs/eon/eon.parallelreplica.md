# {py:mod}`eon.parallelreplica`

```{py:module} eon.parallelreplica
```

```{autodoc2-docstring} eon.parallelreplica
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`parallelreplica <eon.parallelreplica.parallelreplica>`
  - ```{autodoc2-docstring} eon.parallelreplica.parallelreplica
    :summary:
    ```
* - {py:obj}`step <eon.parallelreplica.step>`
  - ```{autodoc2-docstring} eon.parallelreplica.step
    :summary:
    ```
* - {py:obj}`get_statelist <eon.parallelreplica.get_statelist>`
  - ```{autodoc2-docstring} eon.parallelreplica.get_statelist
    :summary:
    ```
* - {py:obj}`get_pr_metadata <eon.parallelreplica.get_pr_metadata>`
  - ```{autodoc2-docstring} eon.parallelreplica.get_pr_metadata
    :summary:
    ```
* - {py:obj}`write_pr_metadata <eon.parallelreplica.write_pr_metadata>`
  - ```{autodoc2-docstring} eon.parallelreplica.write_pr_metadata
    :summary:
    ```
* - {py:obj}`make_searches <eon.parallelreplica.make_searches>`
  - ```{autodoc2-docstring} eon.parallelreplica.make_searches
    :summary:
    ```
* - {py:obj}`register_results <eon.parallelreplica.register_results>`
  - ```{autodoc2-docstring} eon.parallelreplica.register_results
    :summary:
    ```
* - {py:obj}`main <eon.parallelreplica.main>`
  - ```{autodoc2-docstring} eon.parallelreplica.main
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.parallelreplica.logger>`
  - ```{autodoc2-docstring} eon.parallelreplica.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.parallelreplica.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.parallelreplica.logger
```

````

````{py:function} parallelreplica()
:canonical: eon.parallelreplica.parallelreplica

```{autodoc2-docstring} eon.parallelreplica.parallelreplica
```
````

````{py:function} step(current_time, current_state, states, transition)
:canonical: eon.parallelreplica.step

```{autodoc2-docstring} eon.parallelreplica.step
```
````

````{py:function} get_statelist()
:canonical: eon.parallelreplica.get_statelist

```{autodoc2-docstring} eon.parallelreplica.get_statelist
```
````

````{py:function} get_pr_metadata()
:canonical: eon.parallelreplica.get_pr_metadata

```{autodoc2-docstring} eon.parallelreplica.get_pr_metadata
```
````

````{py:function} write_pr_metadata(parser, current_state_num, time, wuid)
:canonical: eon.parallelreplica.write_pr_metadata

```{autodoc2-docstring} eon.parallelreplica.write_pr_metadata
```
````

````{py:function} make_searches(comm, current_state, wuid)
:canonical: eon.parallelreplica.make_searches

```{autodoc2-docstring} eon.parallelreplica.make_searches
```
````

````{py:function} register_results(comm, current_state, states)
:canonical: eon.parallelreplica.register_results

```{autodoc2-docstring} eon.parallelreplica.register_results
```
````

````{py:function} main()
:canonical: eon.parallelreplica.main

```{autodoc2-docstring} eon.parallelreplica.main
```
````
