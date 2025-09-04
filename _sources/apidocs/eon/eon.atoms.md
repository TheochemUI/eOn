# {py:mod}`eon.atoms`

```{py:module} eon.atoms
```

```{autodoc2-docstring} eon.atoms
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Atoms <eon.atoms.Atoms>`
  - ```{autodoc2-docstring} eon.atoms.Atoms
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`pbc <eon.atoms.pbc>`
  - ```{autodoc2-docstring} eon.atoms.pbc
    :summary:
    ```
* - {py:obj}`per_atom_norm <eon.atoms.per_atom_norm>`
  - ```{autodoc2-docstring} eon.atoms.per_atom_norm
    :summary:
    ```
* - {py:obj}`per_atom_norm_gen <eon.atoms.per_atom_norm_gen>`
  - ```{autodoc2-docstring} eon.atoms.per_atom_norm_gen
    :summary:
    ```
* - {py:obj}`get_process_atoms <eon.atoms.get_process_atoms>`
  - ```{autodoc2-docstring} eon.atoms.get_process_atoms
    :summary:
    ```
* - {py:obj}`identical <eon.atoms.identical>`
  - ```{autodoc2-docstring} eon.atoms.identical
    :summary:
    ```
* - {py:obj}`brute_neighbor_list <eon.atoms.brute_neighbor_list>`
  - ```{autodoc2-docstring} eon.atoms.brute_neighbor_list
    :summary:
    ```
* - {py:obj}`sweep_and_prune <eon.atoms.sweep_and_prune>`
  - ```{autodoc2-docstring} eon.atoms.sweep_and_prune
    :summary:
    ```
* - {py:obj}`neighbor_list <eon.atoms.neighbor_list>`
  - ```{autodoc2-docstring} eon.atoms.neighbor_list
    :summary:
    ```
* - {py:obj}`neighbor_list_vectors <eon.atoms.neighbor_list_vectors>`
  - ```{autodoc2-docstring} eon.atoms.neighbor_list_vectors
    :summary:
    ```
* - {py:obj}`coordination_numbers <eon.atoms.coordination_numbers>`
  - ```{autodoc2-docstring} eon.atoms.coordination_numbers
    :summary:
    ```
* - {py:obj}`least_coordinated <eon.atoms.least_coordinated>`
  - ```{autodoc2-docstring} eon.atoms.least_coordinated
    :summary:
    ```
* - {py:obj}`match <eon.atoms.match>`
  - ```{autodoc2-docstring} eon.atoms.match
    :summary:
    ```
* - {py:obj}`point_energy_match <eon.atoms.point_energy_match>`
  - ```{autodoc2-docstring} eon.atoms.point_energy_match
    :summary:
    ```
* - {py:obj}`points_energies_match <eon.atoms.points_energies_match>`
  - ```{autodoc2-docstring} eon.atoms.points_energies_match
    :summary:
    ```
* - {py:obj}`rot_match <eon.atoms.rot_match>`
  - ```{autodoc2-docstring} eon.atoms.rot_match
    :summary:
    ```
* - {py:obj}`rotm <eon.atoms.rotm>`
  - ```{autodoc2-docstring} eon.atoms.rotm
    :summary:
    ```
* - {py:obj}`cna <eon.atoms.cna>`
  - ```{autodoc2-docstring} eon.atoms.cna
    :summary:
    ```
* - {py:obj}`not_HCP_or_FCC <eon.atoms.not_HCP_or_FCC>`
  - ```{autodoc2-docstring} eon.atoms.not_HCP_or_FCC
    :summary:
    ```
* - {py:obj}`cnat <eon.atoms.cnat>`
  - ```{autodoc2-docstring} eon.atoms.cnat
    :summary:
    ```
* - {py:obj}`cnar <eon.atoms.cnar>`
  - ```{autodoc2-docstring} eon.atoms.cnar
    :summary:
    ```
* - {py:obj}`not_TCP <eon.atoms.not_TCP>`
  - ```{autodoc2-docstring} eon.atoms.not_TCP
    :summary:
    ```
* - {py:obj}`not_TCP_or_BCC <eon.atoms.not_TCP_or_BCC>`
  - ```{autodoc2-docstring} eon.atoms.not_TCP_or_BCC
    :summary:
    ```
* - {py:obj}`get_mappings <eon.atoms.get_mappings>`
  - ```{autodoc2-docstring} eon.atoms.get_mappings
    :summary:
    ```
* - {py:obj}`get_rotation_matrix <eon.atoms.get_rotation_matrix>`
  - ```{autodoc2-docstring} eon.atoms.get_rotation_matrix
    :summary:
    ```
* - {py:obj}`rotate <eon.atoms.rotate>`
  - ```{autodoc2-docstring} eon.atoms.rotate
    :summary:
    ```
* - {py:obj}`internal_motion <eon.atoms.internal_motion>`
  - ```{autodoc2-docstring} eon.atoms.internal_motion
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.atoms.logger>`
  - ```{autodoc2-docstring} eon.atoms.logger
    :summary:
    ```
* - {py:obj}`elements <eon.atoms.elements>`
  - ```{autodoc2-docstring} eon.atoms.elements
    :summary:
    ```
* - {py:obj}`numElements <eon.atoms.numElements>`
  - ```{autodoc2-docstring} eon.atoms.numElements
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.atoms.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.atoms.logger
```

````

`````{py:class} Atoms(n_atoms)
:canonical: eon.atoms.Atoms

```{autodoc2-docstring} eon.atoms.Atoms
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.atoms.Atoms.__init__
```

````{py:method} copy()
:canonical: eon.atoms.Atoms.copy

```{autodoc2-docstring} eon.atoms.Atoms.copy
```

````

````{py:method} free_r()
:canonical: eon.atoms.Atoms.free_r

```{autodoc2-docstring} eon.atoms.Atoms.free_r
```

````

````{py:method} append(r, free, name, mass)
:canonical: eon.atoms.Atoms.append

```{autodoc2-docstring} eon.atoms.Atoms.append
```

````

`````

````{py:function} pbc(r, box, ibox=None)
:canonical: eon.atoms.pbc

```{autodoc2-docstring} eon.atoms.pbc
```
````

````{py:function} per_atom_norm(v, box, ibox=None)
:canonical: eon.atoms.per_atom_norm

```{autodoc2-docstring} eon.atoms.per_atom_norm
```
````

````{py:function} per_atom_norm_gen(v, box, ibox=None)
:canonical: eon.atoms.per_atom_norm_gen

```{autodoc2-docstring} eon.atoms.per_atom_norm_gen
```
````

````{py:function} get_process_atoms(r, p, epsilon_r=0.2, nshells=1)
:canonical: eon.atoms.get_process_atoms

```{autodoc2-docstring} eon.atoms.get_process_atoms
```
````

````{py:function} identical(atoms1, atoms2)
:canonical: eon.atoms.identical

```{autodoc2-docstring} eon.atoms.identical
```
````

````{py:function} brute_neighbor_list(p, cutoff)
:canonical: eon.atoms.brute_neighbor_list

```{autodoc2-docstring} eon.atoms.brute_neighbor_list
```
````

````{py:function} sweep_and_prune(p_in, cutoff, strict=True, bc=True)
:canonical: eon.atoms.sweep_and_prune

```{autodoc2-docstring} eon.atoms.sweep_and_prune
```
````

````{py:function} neighbor_list(p, cutoff, brute=False)
:canonical: eon.atoms.neighbor_list

```{autodoc2-docstring} eon.atoms.neighbor_list
```
````

````{py:function} neighbor_list_vectors(p, cutoff, brute=False)
:canonical: eon.atoms.neighbor_list_vectors

```{autodoc2-docstring} eon.atoms.neighbor_list_vectors
```
````

````{py:function} coordination_numbers(p, cutoff, brute=False)
:canonical: eon.atoms.coordination_numbers

```{autodoc2-docstring} eon.atoms.coordination_numbers
```
````

````{py:function} least_coordinated(p, cutoff, brute=False)
:canonical: eon.atoms.least_coordinated

```{autodoc2-docstring} eon.atoms.least_coordinated
```
````

````{py:function} match(a, b, eps_r, neighbor_cutoff, indistinguishable)
:canonical: eon.atoms.match

```{autodoc2-docstring} eon.atoms.match
```
````

````{py:function} point_energy_match(file_a, energy_a, file_b, energy_b)
:canonical: eon.atoms.point_energy_match

```{autodoc2-docstring} eon.atoms.point_energy_match
```
````

````{py:function} points_energies_match(file_a, energy_a, files_b, energies_b)
:canonical: eon.atoms.points_energies_match

```{autodoc2-docstring} eon.atoms.points_energies_match
```
````

````{py:function} rot_match(a, b)
:canonical: eon.atoms.rot_match

```{autodoc2-docstring} eon.atoms.rot_match
```
````

````{py:function} rotm(axis, theta)
:canonical: eon.atoms.rotm

```{autodoc2-docstring} eon.atoms.rotm
```
````

````{py:function} cna(p, cutoff, brute=False)
:canonical: eon.atoms.cna

```{autodoc2-docstring} eon.atoms.cna
```
````

````{py:function} not_HCP_or_FCC(p, cutoff, brute=False)
:canonical: eon.atoms.not_HCP_or_FCC

```{autodoc2-docstring} eon.atoms.not_HCP_or_FCC
```
````

````{py:function} cnat(p, cutoff, brute=False)
:canonical: eon.atoms.cnat

```{autodoc2-docstring} eon.atoms.cnat
```
````

````{py:function} cnar(p, cutoff, brute=False)
:canonical: eon.atoms.cnar

```{autodoc2-docstring} eon.atoms.cnar
```
````

````{py:function} not_TCP(p, cutoff, brute=False)
:canonical: eon.atoms.not_TCP

```{autodoc2-docstring} eon.atoms.not_TCP
```
````

````{py:function} not_TCP_or_BCC(p, cutoff, brute=False)
:canonical: eon.atoms.not_TCP_or_BCC

```{autodoc2-docstring} eon.atoms.not_TCP_or_BCC
```
````

````{py:function} get_mappings(a, b, eps_r, neighbor_cutoff, mappings=None)
:canonical: eon.atoms.get_mappings

```{autodoc2-docstring} eon.atoms.get_mappings
```
````

````{py:function} get_rotation_matrix(axis, theta)
:canonical: eon.atoms.get_rotation_matrix

```{autodoc2-docstring} eon.atoms.get_rotation_matrix
```
````

````{py:function} rotate(r, axis, center, angle)
:canonical: eon.atoms.rotate

```{autodoc2-docstring} eon.atoms.rotate
```
````

````{py:function} internal_motion(a, b)
:canonical: eon.atoms.internal_motion

```{autodoc2-docstring} eon.atoms.internal_motion
```
````

````{py:data} elements
:canonical: eon.atoms.elements
:value: >
   None

```{autodoc2-docstring} eon.atoms.elements
```

````

````{py:data} numElements
:canonical: eon.atoms.numElements
:value: >
   119

```{autodoc2-docstring} eon.atoms.numElements
```

````
