# {py:mod}`eon.fileio`

```{py:module} eon.fileio
```

```{autodoc2-docstring} eon.fileio
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`ini <eon.fileio.ini>`
  -
* - {py:obj}`Dynamics <eon.fileio.Dynamics>`
  - ```{autodoc2-docstring} eon.fileio.Dynamics
    :summary:
    ```
* - {py:obj}`Table <eon.fileio.Table>`
  - ```{autodoc2-docstring} eon.fileio.Table
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`save_prng_state <eon.fileio.save_prng_state>`
  - ```{autodoc2-docstring} eon.fileio.save_prng_state
    :summary:
    ```
* - {py:obj}`get_prng_state <eon.fileio.get_prng_state>`
  - ```{autodoc2-docstring} eon.fileio.get_prng_state
    :summary:
    ```
* - {py:obj}`length_angle_to_box <eon.fileio.length_angle_to_box>`
  - ```{autodoc2-docstring} eon.fileio.length_angle_to_box
    :summary:
    ```
* - {py:obj}`box_to_length_angle <eon.fileio.box_to_length_angle>`
  - ```{autodoc2-docstring} eon.fileio.box_to_length_angle
    :summary:
    ```
* - {py:obj}`loadcons <eon.fileio.loadcons>`
  - ```{autodoc2-docstring} eon.fileio.loadcons
    :summary:
    ```
* - {py:obj}`loadposcars <eon.fileio.loadposcars>`
  - ```{autodoc2-docstring} eon.fileio.loadposcars
    :summary:
    ```
* - {py:obj}`loadcon <eon.fileio.loadcon>`
  - ```{autodoc2-docstring} eon.fileio.loadcon
    :summary:
    ```
* - {py:obj}`savecon <eon.fileio.savecon>`
  - ```{autodoc2-docstring} eon.fileio.savecon
    :summary:
    ```
* - {py:obj}`load_mode <eon.fileio.load_mode>`
  - ```{autodoc2-docstring} eon.fileio.load_mode
    :summary:
    ```
* - {py:obj}`save_mode <eon.fileio.save_mode>`
  - ```{autodoc2-docstring} eon.fileio.save_mode
    :summary:
    ```
* - {py:obj}`save_results_dat <eon.fileio.save_results_dat>`
  - ```{autodoc2-docstring} eon.fileio.save_results_dat
    :summary:
    ```
* - {py:obj}`modify_config <eon.fileio.modify_config>`
  - ```{autodoc2-docstring} eon.fileio.modify_config
    :summary:
    ```
* - {py:obj}`parse_results <eon.fileio.parse_results>`
  - ```{autodoc2-docstring} eon.fileio.parse_results
    :summary:
    ```
* - {py:obj}`loadposcar <eon.fileio.loadposcar>`
  - ```{autodoc2-docstring} eon.fileio.loadposcar
    :summary:
    ```
* - {py:obj}`saveposcar <eon.fileio.saveposcar>`
  - ```{autodoc2-docstring} eon.fileio.saveposcar
    :summary:
    ```
* - {py:obj}`load_potfiles <eon.fileio.load_potfiles>`
  - ```{autodoc2-docstring} eon.fileio.load_potfiles
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.fileio.logger>`
  - ```{autodoc2-docstring} eon.fileio.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.fileio.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.fileio.logger
```

````

````{py:function} save_prng_state()
:canonical: eon.fileio.save_prng_state

```{autodoc2-docstring} eon.fileio.save_prng_state
```
````

````{py:function} get_prng_state()
:canonical: eon.fileio.get_prng_state

```{autodoc2-docstring} eon.fileio.get_prng_state
```
````

````{py:function} length_angle_to_box(boxlengths, angles)
:canonical: eon.fileio.length_angle_to_box

```{autodoc2-docstring} eon.fileio.length_angle_to_box
```
````

````{py:function} box_to_length_angle(box)
:canonical: eon.fileio.box_to_length_angle

```{autodoc2-docstring} eon.fileio.box_to_length_angle
```
````

````{py:function} loadcons(filename)
:canonical: eon.fileio.loadcons

```{autodoc2-docstring} eon.fileio.loadcons
```
````

````{py:function} loadposcars(filename)
:canonical: eon.fileio.loadposcars

```{autodoc2-docstring} eon.fileio.loadposcars
```
````

````{py:function} loadcon(filein, reset=True)
:canonical: eon.fileio.loadcon

```{autodoc2-docstring} eon.fileio.loadcon
```
````

````{py:function} savecon(fileout, p, w='w')
:canonical: eon.fileio.savecon

```{autodoc2-docstring} eon.fileio.savecon
```
````

````{py:function} load_mode(modefilein)
:canonical: eon.fileio.load_mode

```{autodoc2-docstring} eon.fileio.load_mode
```
````

````{py:function} save_mode(modefileout, displace_vector)
:canonical: eon.fileio.save_mode

```{autodoc2-docstring} eon.fileio.save_mode
```
````

````{py:function} save_results_dat(fileout, results)
:canonical: eon.fileio.save_results_dat

```{autodoc2-docstring} eon.fileio.save_results_dat
```
````

````{py:function} modify_config(config_path, changes)
:canonical: eon.fileio.modify_config

```{autodoc2-docstring} eon.fileio.modify_config
```
````

````{py:function} parse_results(filein)
:canonical: eon.fileio.parse_results

```{autodoc2-docstring} eon.fileio.parse_results
```
````

````{py:function} loadposcar(filein)
:canonical: eon.fileio.loadposcar

```{autodoc2-docstring} eon.fileio.loadposcar
```
````

````{py:function} saveposcar(fileout, p, w='w', direct=False)
:canonical: eon.fileio.saveposcar

```{autodoc2-docstring} eon.fileio.saveposcar
```
````

`````{py:class} ini(filenames)
:canonical: eon.fileio.ini

Bases: {py:obj}`configparser.ConfigParser`

````{py:method} read()
:canonical: eon.fileio.ini.read

````

````{py:method} get(section, option, default='ini_no_default', **kwargs)
:canonical: eon.fileio.ini.get

````

````{py:method} getint(*args)
:canonical: eon.fileio.ini.getint
:abstractmethod:

```{autodoc2-docstring} eon.fileio.ini.getint
```

````

````{py:method} getfloat(*args)
:canonical: eon.fileio.ini.getfloat
:abstractmethod:

```{autodoc2-docstring} eon.fileio.ini.getfloat
```

````

````{py:method} getboolean(*args)
:canonical: eon.fileio.ini.getboolean
:abstractmethod:

```{autodoc2-docstring} eon.fileio.ini.getboolean
```

````

````{py:method} set(section, option, value)
:canonical: eon.fileio.ini.set

````

`````

`````{py:class} Dynamics(filename)
:canonical: eon.fileio.Dynamics

```{autodoc2-docstring} eon.fileio.Dynamics
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.fileio.Dynamics.__init__
```

````{py:method} append(reactant_id, process_id, product_id, step_time, total_time, barrier, rate, energy)
:canonical: eon.fileio.Dynamics.append

```{autodoc2-docstring} eon.fileio.Dynamics.append
```

````

````{py:method} append_sb(reactant_id, process_id, product_id, step_time, total_time, basin_id, rate, energy)
:canonical: eon.fileio.Dynamics.append_sb

```{autodoc2-docstring} eon.fileio.Dynamics.append_sb
```

````

````{py:method} get()
:canonical: eon.fileio.Dynamics.get

```{autodoc2-docstring} eon.fileio.Dynamics.get
```

````

`````

````{py:function} load_potfiles(pot_dir)
:canonical: eon.fileio.load_potfiles

```{autodoc2-docstring} eon.fileio.load_potfiles
```
````

```{py:exception} TableException()
:canonical: eon.fileio.TableException

Bases: {py:obj}`Exception`

```

`````{py:class} Table(filename, columns=None, overwrite=False)
:canonical: eon.fileio.Table

```{autodoc2-docstring} eon.fileio.Table
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.fileio.Table.__init__
```

````{py:method} init()
:canonical: eon.fileio.Table.init

```{autodoc2-docstring} eon.fileio.Table.init
```

````

````{py:method} read(filename)
:canonical: eon.fileio.Table.read

```{autodoc2-docstring} eon.fileio.Table.read
```

````

````{py:method} write()
:canonical: eon.fileio.Table.write

```{autodoc2-docstring} eon.fileio.Table.write
```

````

````{py:method} writefilehandle(filehandle)
:canonical: eon.fileio.Table.writefilehandle

```{autodoc2-docstring} eon.fileio.Table.writefilehandle
```

````

````{py:method} add_row(row)
:canonical: eon.fileio.Table.add_row

```{autodoc2-docstring} eon.fileio.Table.add_row
```

````

````{py:method} delete_row(column, value)
:canonical: eon.fileio.Table.delete_row

```{autodoc2-docstring} eon.fileio.Table.delete_row
```

````

````{py:method} delete_row_func(column, func)
:canonical: eon.fileio.Table.delete_row_func

```{autodoc2-docstring} eon.fileio.Table.delete_row_func
```

````

````{py:method} find_value(column, func)
:canonical: eon.fileio.Table.find_value

```{autodoc2-docstring} eon.fileio.Table.find_value
```

````

````{py:method} find_row(column, func)
:canonical: eon.fileio.Table.find_row

```{autodoc2-docstring} eon.fileio.Table.find_row
```

````

````{py:method} min_value(column)
:canonical: eon.fileio.Table.min_value

```{autodoc2-docstring} eon.fileio.Table.min_value
```

````

````{py:method} min_row(column)
:canonical: eon.fileio.Table.min_row

```{autodoc2-docstring} eon.fileio.Table.min_row
```

````

````{py:method} max_value(column)
:canonical: eon.fileio.Table.max_value

```{autodoc2-docstring} eon.fileio.Table.max_value
```

````

````{py:method} max_row(column)
:canonical: eon.fileio.Table.max_row

```{autodoc2-docstring} eon.fileio.Table.max_row
```

````

````{py:method} get_row(column, value)
:canonical: eon.fileio.Table.get_row

```{autodoc2-docstring} eon.fileio.Table.get_row
```

````

````{py:method} get_rows(column, value)
:canonical: eon.fileio.Table.get_rows

```{autodoc2-docstring} eon.fileio.Table.get_rows
```

````

````{py:method} get_column(column)
:canonical: eon.fileio.Table.get_column

```{autodoc2-docstring} eon.fileio.Table.get_column
```

````

`````
