# {py:mod}`eon.mcamc.mcamc`

```{py:module} eon.mcamc.mcamc
```

```{autodoc2-docstring} eon.mcamc.mcamc
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`estimate_condition <eon.mcamc.mcamc.estimate_condition>`
  - ```{autodoc2-docstring} eon.mcamc.mcamc.estimate_condition
    :summary:
    ```
* - {py:obj}`guess_precision <eon.mcamc.mcamc.guess_precision>`
  - ```{autodoc2-docstring} eon.mcamc.mcamc.guess_precision
    :summary:
    ```
* - {py:obj}`c_mcamc <eon.mcamc.mcamc.c_mcamc>`
  - ```{autodoc2-docstring} eon.mcamc.mcamc.c_mcamc
    :summary:
    ```
* - {py:obj}`np_mcamc <eon.mcamc.mcamc.np_mcamc>`
  - ```{autodoc2-docstring} eon.mcamc.mcamc.np_mcamc
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.mcamc.mcamc.logger>`
  - ```{autodoc2-docstring} eon.mcamc.mcamc.logger
    :summary:
    ```
* - {py:obj}`libpath <eon.mcamc.mcamc.libpath>`
  - ```{autodoc2-docstring} eon.mcamc.mcamc.libpath
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.mcamc.mcamc.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.mcamc.mcamc.logger
```

````

````{py:function} estimate_condition(Q, R)
:canonical: eon.mcamc.mcamc.estimate_condition

```{autodoc2-docstring} eon.mcamc.mcamc.estimate_condition
```
````

````{py:function} guess_precision(Q, R)
:canonical: eon.mcamc.mcamc.guess_precision

```{autodoc2-docstring} eon.mcamc.mcamc.guess_precision
```
````

````{py:function} c_mcamc(Q, R, c, prec='dd')
:canonical: eon.mcamc.mcamc.c_mcamc

```{autodoc2-docstring} eon.mcamc.mcamc.c_mcamc
```
````

````{py:function} np_mcamc(Q, R, c, prec='NA')
:canonical: eon.mcamc.mcamc.np_mcamc

```{autodoc2-docstring} eon.mcamc.mcamc.np_mcamc
```
````

````{py:data} libpath
:canonical: eon.mcamc.mcamc.libpath
:value: >
   'join(...)'

```{autodoc2-docstring} eon.mcamc.mcamc.libpath
```

````
