# {py:mod}`eon.akmc`

```{py:module} eon.akmc
```

```{autodoc2-docstring} eon.akmc
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`akmc <eon.akmc.akmc>`
  - ```{autodoc2-docstring} eon.akmc.akmc
    :summary:
    ```
* - {py:obj}`get_akmc_metadata <eon.akmc.get_akmc_metadata>`
  - ```{autodoc2-docstring} eon.akmc.get_akmc_metadata
    :summary:
    ```
* - {py:obj}`write_akmc_metadata <eon.akmc.write_akmc_metadata>`
  - ```{autodoc2-docstring} eon.akmc.write_akmc_metadata
    :summary:
    ```
* - {py:obj}`get_statelist <eon.akmc.get_statelist>`
  - ```{autodoc2-docstring} eon.akmc.get_statelist
    :summary:
    ```
* - {py:obj}`get_superbasin_scheme <eon.akmc.get_superbasin_scheme>`
  - ```{autodoc2-docstring} eon.akmc.get_superbasin_scheme
    :summary:
    ```
* - {py:obj}`kmc_step <eon.akmc.kmc_step>`
  - ```{autodoc2-docstring} eon.akmc.kmc_step
    :summary:
    ```
* - {py:obj}`main <eon.akmc.main>`
  - ```{autodoc2-docstring} eon.akmc.main
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.akmc.logger>`
  - ```{autodoc2-docstring} eon.akmc.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.akmc.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.akmc.logger
```

````

````{py:function} akmc(config: eon.config.ConfigClass = EON_CONFIG, steps=0)
:canonical: eon.akmc.akmc

```{autodoc2-docstring} eon.akmc.akmc
```
````

````{py:function} get_akmc_metadata(config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.akmc.get_akmc_metadata

```{autodoc2-docstring} eon.akmc.get_akmc_metadata
```
````

````{py:function} write_akmc_metadata(parser, current_state_num, time, previous_state_num, previous_temperature)
:canonical: eon.akmc.write_akmc_metadata

```{autodoc2-docstring} eon.akmc.write_akmc_metadata
```
````

````{py:function} get_statelist(kT, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.akmc.get_statelist

```{autodoc2-docstring} eon.akmc.get_statelist
```
````

````{py:function} get_superbasin_scheme(states, config)
:canonical: eon.akmc.get_superbasin_scheme

```{autodoc2-docstring} eon.akmc.get_superbasin_scheme
```
````

````{py:function} kmc_step(current_state, states, time, kT, superbasining, steps=0, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.akmc.kmc_step

```{autodoc2-docstring} eon.akmc.kmc_step
```
````

````{py:function} main(config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.akmc.main

```{autodoc2-docstring} eon.akmc.main
```
````
