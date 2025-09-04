# {py:mod}`eon.movie`

```{py:module} eon.movie
```

```{autodoc2-docstring} eon.movie
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Graph <eon.movie.Graph>`
  - ```{autodoc2-docstring} eon.movie.Graph
    :summary:
    ```
* - {py:obj}`priorityDictionary <eon.movie.priorityDictionary>`
  -
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`make_movie <eon.movie.make_movie>`
  - ```{autodoc2-docstring} eon.movie.make_movie
    :summary:
    ```
* - {py:obj}`get_trajectory <eon.movie.get_trajectory>`
  - ```{autodoc2-docstring} eon.movie.get_trajectory
    :summary:
    ```
* - {py:obj}`processes <eon.movie.processes>`
  - ```{autodoc2-docstring} eon.movie.processes
    :summary:
    ```
* - {py:obj}`dot <eon.movie.dot>`
  - ```{autodoc2-docstring} eon.movie.dot
    :summary:
    ```
* - {py:obj}`dynamics <eon.movie.dynamics>`
  - ```{autodoc2-docstring} eon.movie.dynamics
    :summary:
    ```
* - {py:obj}`make_graph <eon.movie.make_graph>`
  - ```{autodoc2-docstring} eon.movie.make_graph
    :summary:
    ```
* - {py:obj}`fastest_path <eon.movie.fastest_path>`
  - ```{autodoc2-docstring} eon.movie.fastest_path
    :summary:
    ```
* - {py:obj}`get_fastest_process_id <eon.movie.get_fastest_process_id>`
  - ```{autodoc2-docstring} eon.movie.get_fastest_process_id
    :summary:
    ```
* - {py:obj}`get_fastest_process_rate <eon.movie.get_fastest_process_rate>`
  - ```{autodoc2-docstring} eon.movie.get_fastest_process_rate
    :summary:
    ```
````

### API

````{py:function} make_movie(movie_type, path_root, states, separate_files=False)
:canonical: eon.movie.make_movie

```{autodoc2-docstring} eon.movie.make_movie
```
````

````{py:function} get_trajectory(trajectory_path)
:canonical: eon.movie.get_trajectory

```{autodoc2-docstring} eon.movie.get_trajectory
```
````

````{py:function} processes(states, statenr, limit)
:canonical: eon.movie.processes

```{autodoc2-docstring} eon.movie.processes
```
````

````{py:function} dot(path_root, states)
:canonical: eon.movie.dot

```{autodoc2-docstring} eon.movie.dot
```
````

````{py:function} dynamics(path_root, states, unique=False)
:canonical: eon.movie.dynamics

```{autodoc2-docstring} eon.movie.dynamics
```
````

````{py:function} make_graph(states)
:canonical: eon.movie.make_graph

```{autodoc2-docstring} eon.movie.make_graph
```
````

````{py:function} fastest_path(path_root, states, full=False)
:canonical: eon.movie.fastest_path

```{autodoc2-docstring} eon.movie.fastest_path
```
````

````{py:function} get_fastest_process_id(state1, state2)
:canonical: eon.movie.get_fastest_process_id

```{autodoc2-docstring} eon.movie.get_fastest_process_id
```
````

````{py:function} get_fastest_process_rate(state1, state2)
:canonical: eon.movie.get_fastest_process_rate

```{autodoc2-docstring} eon.movie.get_fastest_process_rate
```
````

`````{py:class} Graph(name='')
:canonical: eon.movie.Graph

```{autodoc2-docstring} eon.movie.Graph
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.movie.Graph.__init__
```

````{py:method} dot()
:canonical: eon.movie.Graph.dot

```{autodoc2-docstring} eon.movie.Graph.dot
```

````

````{py:method} add_node(node)
:canonical: eon.movie.Graph.add_node

```{autodoc2-docstring} eon.movie.Graph.add_node
```

````

````{py:method} add_edge(node1, node2, weight=1.0)
:canonical: eon.movie.Graph.add_edge

```{autodoc2-docstring} eon.movie.Graph.add_edge
```

````

````{py:method} neighbors(node)
:canonical: eon.movie.Graph.neighbors

```{autodoc2-docstring} eon.movie.Graph.neighbors
```

````

````{py:method} nodes()
:canonical: eon.movie.Graph.nodes

```{autodoc2-docstring} eon.movie.Graph.nodes
```

````

````{py:method} dijkstra(start, end=None)
:canonical: eon.movie.Graph.dijkstra

```{autodoc2-docstring} eon.movie.Graph.dijkstra
```

````

````{py:method} shortest_path(start, end)
:canonical: eon.movie.Graph.shortest_path

```{autodoc2-docstring} eon.movie.Graph.shortest_path
```

````

`````

`````{py:class} priorityDictionary()
:canonical: eon.movie.priorityDictionary

Bases: {py:obj}`dict`

````{py:method} smallest()
:canonical: eon.movie.priorityDictionary.smallest

```{autodoc2-docstring} eon.movie.priorityDictionary.smallest
```

````

````{py:method} setdefault(key, val)
:canonical: eon.movie.priorityDictionary.setdefault

```{autodoc2-docstring} eon.movie.priorityDictionary.setdefault
```

````

`````
