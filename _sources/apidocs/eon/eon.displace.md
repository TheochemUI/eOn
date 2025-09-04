# {py:mod}`eon.displace`

```{py:module} eon.displace
```

```{autodoc2-docstring} eon.displace
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`DisplacementManager <eon.displace.DisplacementManager>`
  - ```{autodoc2-docstring} eon.displace.DisplacementManager
    :summary:
    ```
* - {py:obj}`Displace <eon.displace.Displace>`
  - ```{autodoc2-docstring} eon.displace.Displace
    :summary:
    ```
* - {py:obj}`Undercoordinated <eon.displace.Undercoordinated>`
  - ```{autodoc2-docstring} eon.displace.Undercoordinated
    :summary:
    ```
* - {py:obj}`Leastcoordinated <eon.displace.Leastcoordinated>`
  - ```{autodoc2-docstring} eon.displace.Leastcoordinated
    :summary:
    ```
* - {py:obj}`ListedAtoms <eon.displace.ListedAtoms>`
  - ```{autodoc2-docstring} eon.displace.ListedAtoms
    :summary:
    ```
* - {py:obj}`ListedTypes <eon.displace.ListedTypes>`
  - ```{autodoc2-docstring} eon.displace.ListedTypes
    :summary:
    ```
* - {py:obj}`Random <eon.displace.Random>`
  - ```{autodoc2-docstring} eon.displace.Random
    :summary:
    ```
* - {py:obj}`NotFCCorHCP <eon.displace.NotFCCorHCP>`
  - ```{autodoc2-docstring} eon.displace.NotFCCorHCP
    :summary:
    ```
* - {py:obj}`NotTCPorBCC <eon.displace.NotTCPorBCC>`
  - ```{autodoc2-docstring} eon.displace.NotTCPorBCC
    :summary:
    ```
* - {py:obj}`NotTCP <eon.displace.NotTCP>`
  - ```{autodoc2-docstring} eon.displace.NotTCP
    :summary:
    ```
* - {py:obj}`Water <eon.displace.Water>`
  - ```{autodoc2-docstring} eon.displace.Water
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.displace.logger>`
  - ```{autodoc2-docstring} eon.displace.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.displace.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.displace.logger
```

````

`````{py:class} DisplacementManager(reactant, moved_atoms, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.DisplacementManager

```{autodoc2-docstring} eon.displace.DisplacementManager
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.DisplacementManager.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.DisplacementManager.make_displacement

```{autodoc2-docstring} eon.displace.DisplacementManager.make_displacement
```

````

`````

```{py:exception} NotImplementedError()
:canonical: eon.displace.NotImplementedError

Bases: {py:obj}`Exception`

```

```{py:exception} DisplaceError()
:canonical: eon.displace.DisplaceError

Bases: {py:obj}`Exception`

```

`````{py:class} Displace(reactant, std_dev, radius, hole_epicenters, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.Displace

```{autodoc2-docstring} eon.displace.Displace
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.Displace.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.Displace.make_displacement
:abstractmethod:

```{autodoc2-docstring} eon.displace.Displace.make_displacement
```

````

````{py:method} get_displacement(atom_index)
:canonical: eon.displace.Displace.get_displacement

```{autodoc2-docstring} eon.displace.Displace.get_displacement
```

````

````{py:method} filter_epicenters(epicenters)
:canonical: eon.displace.Displace.filter_epicenters

```{autodoc2-docstring} eon.displace.Displace.filter_epicenters
```

````

`````

`````{py:class} Undercoordinated(reactant, max_coordination, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.Undercoordinated

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.Undercoordinated
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.Undercoordinated.__init__
```

````{py:method} init()
:canonical: eon.displace.Undercoordinated.init

```{autodoc2-docstring} eon.displace.Undercoordinated.init
```

````

````{py:method} make_displacement()
:canonical: eon.displace.Undercoordinated.make_displacement

```{autodoc2-docstring} eon.displace.Undercoordinated.make_displacement
```

````

`````

`````{py:class} Leastcoordinated(reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.Leastcoordinated

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.Leastcoordinated
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.Leastcoordinated.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.Leastcoordinated.make_displacement

```{autodoc2-docstring} eon.displace.Leastcoordinated.make_displacement
```

````

`````

`````{py:class} ListedAtoms(reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, displace_all=False, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.ListedAtoms

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.ListedAtoms
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.ListedAtoms.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.ListedAtoms.make_displacement

```{autodoc2-docstring} eon.displace.ListedAtoms.make_displacement
```

````

`````

`````{py:class} ListedTypes(reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, displace_all=False, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.ListedTypes

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.ListedTypes
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.ListedTypes.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.ListedTypes.make_displacement

```{autodoc2-docstring} eon.displace.ListedTypes.make_displacement
```

````

`````

`````{py:class} Random(reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.Random

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.Random
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.Random.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.Random.make_displacement

```{autodoc2-docstring} eon.displace.Random.make_displacement
```

````

`````

`````{py:class} NotFCCorHCP(reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.NotFCCorHCP

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.NotFCCorHCP
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.NotFCCorHCP.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.NotFCCorHCP.make_displacement

```{autodoc2-docstring} eon.displace.NotFCCorHCP.make_displacement
```

````

`````

`````{py:class} NotTCPorBCC(reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.NotTCPorBCC

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.NotTCPorBCC
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.NotTCPorBCC.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.NotTCPorBCC.make_displacement

```{autodoc2-docstring} eon.displace.NotTCPorBCC.make_displacement
```

````

`````

`````{py:class} NotTCP(reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.displace.NotTCP

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.NotTCP
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.NotTCP.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.NotTCP.make_displacement

```{autodoc2-docstring} eon.displace.NotTCP.make_displacement
```

````

`````

`````{py:class} Water(reactant, stdev_translation, stdev_rotation, molecule_list=[], random=0)
:canonical: eon.displace.Water

Bases: {py:obj}`eon.displace.Displace`

```{autodoc2-docstring} eon.displace.Water
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.displace.Water.__init__
```

````{py:method} make_displacement()
:canonical: eon.displace.Water.make_displacement

```{autodoc2-docstring} eon.displace.Water.make_displacement
```

````

````{py:method} get_displacement()
:canonical: eon.displace.Water.get_displacement

```{autodoc2-docstring} eon.displace.Water.get_displacement
```

````

````{py:method} rotate_water(hydrogen1, hydrogen2, oxygen, psi, theta, phi, hydrogen_mass=1.0, oxygen_mass=16.0)
:canonical: eon.displace.Water.rotate_water
:staticmethod:

```{autodoc2-docstring} eon.displace.Water.rotate_water
```

````

`````
