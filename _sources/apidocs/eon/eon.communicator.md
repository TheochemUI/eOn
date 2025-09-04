# {py:mod}`eon.communicator`

```{py:module} eon.communicator
```

```{autodoc2-docstring} eon.communicator
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Communicator <eon.communicator.Communicator>`
  - ```{autodoc2-docstring} eon.communicator.Communicator
    :summary:
    ```
* - {py:obj}`MPI <eon.communicator.MPI>`
  - ```{autodoc2-docstring} eon.communicator.MPI
    :summary:
    ```
* - {py:obj}`Local <eon.communicator.Local>`
  - ```{autodoc2-docstring} eon.communicator.Local
    :summary:
    ```
* - {py:obj}`Script <eon.communicator.Script>`
  - ```{autodoc2-docstring} eon.communicator.Script
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`tryint <eon.communicator.tryint>`
  - ```{autodoc2-docstring} eon.communicator.tryint
    :summary:
    ```
* - {py:obj}`alphanum_key <eon.communicator.alphanum_key>`
  - ```{autodoc2-docstring} eon.communicator.alphanum_key
    :summary:
    ```
* - {py:obj}`sort_nicely <eon.communicator.sort_nicely>`
  - ```{autodoc2-docstring} eon.communicator.sort_nicely
    :summary:
    ```
* - {py:obj}`get_communicator <eon.communicator.get_communicator>`
  - ```{autodoc2-docstring} eon.communicator.get_communicator
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`logger <eon.communicator.logger>`
  - ```{autodoc2-docstring} eon.communicator.logger
    :summary:
    ```
````

### API

````{py:data} logger
:canonical: eon.communicator.logger
:value: >
   'getLogger(...)'

```{autodoc2-docstring} eon.communicator.logger
```

````

````{py:function} tryint(s)
:canonical: eon.communicator.tryint

```{autodoc2-docstring} eon.communicator.tryint
```
````

````{py:function} alphanum_key(s)
:canonical: eon.communicator.alphanum_key

```{autodoc2-docstring} eon.communicator.alphanum_key
```
````

````{py:function} sort_nicely(l)
:canonical: eon.communicator.sort_nicely

```{autodoc2-docstring} eon.communicator.sort_nicely
```
````

````{py:function} get_communicator(config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.communicator.get_communicator

```{autodoc2-docstring} eon.communicator.get_communicator
```
````

```{py:exception} NotImplementedError()
:canonical: eon.communicator.NotImplementedError

Bases: {py:obj}`Exception`

```

```{py:exception} CommunicatorError()
:canonical: eon.communicator.CommunicatorError

Bases: {py:obj}`Exception`

```

````{py:exception} EONClientError()
:canonical: eon.communicator.EONClientError

Bases: {py:obj}`Exception`

```{autodoc2-docstring} eon.communicator.EONClientError
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.communicator.EONClientError.__init__
```

````

`````{py:class} Communicator(scratchpath, bundle_size=1, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.communicator.Communicator

```{autodoc2-docstring} eon.communicator.Communicator
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.communicator.Communicator.__init__
```

````{py:method} submit_jobs(data, invariants)
:canonical: eon.communicator.Communicator.submit_jobs
:abstractmethod:

```{autodoc2-docstring} eon.communicator.Communicator.submit_jobs
```

````

````{py:method} get_results(results)
:canonical: eon.communicator.Communicator.get_results
:abstractmethod:

```{autodoc2-docstring} eon.communicator.Communicator.get_results
```

````

````{py:method} get_queue_size()
:canonical: eon.communicator.Communicator.get_queue_size
:abstractmethod:

```{autodoc2-docstring} eon.communicator.Communicator.get_queue_size
```

````

````{py:method} cancel_state(statenumber)
:canonical: eon.communicator.Communicator.cancel_state
:abstractmethod:

```{autodoc2-docstring} eon.communicator.Communicator.cancel_state
```

````

````{py:method} get_bundle_size(job_path)
:canonical: eon.communicator.Communicator.get_bundle_size

```{autodoc2-docstring} eon.communicator.Communicator.get_bundle_size
```

````

````{py:method} unbundle(resultpath, keep_result)
:canonical: eon.communicator.Communicator.unbundle

```{autodoc2-docstring} eon.communicator.Communicator.unbundle
```

````

````{py:method} make_bundles(data, invariants)
:canonical: eon.communicator.Communicator.make_bundles

```{autodoc2-docstring} eon.communicator.Communicator.make_bundles
```

````

`````

`````{py:class} MPI(scratchpath, bundle_size, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.communicator.MPI

Bases: {py:obj}`eon.communicator.Communicator`

```{autodoc2-docstring} eon.communicator.MPI
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.communicator.MPI.__init__
```

````{py:method} submit_jobs(data, invariants)
:canonical: eon.communicator.MPI.submit_jobs

````

````{py:method} run_resume_jobs()
:canonical: eon.communicator.MPI.run_resume_jobs

```{autodoc2-docstring} eon.communicator.MPI.run_resume_jobs
```

````

````{py:method} get_ready_ranks()
:canonical: eon.communicator.MPI.get_ready_ranks

```{autodoc2-docstring} eon.communicator.MPI.get_ready_ranks
```

````

````{py:method} get_queue_size()
:canonical: eon.communicator.MPI.get_queue_size

````

````{py:method} get_results(resultspath, keep_result)
:canonical: eon.communicator.MPI.get_results

````

````{py:method} get_number_in_progress()
:canonical: eon.communicator.MPI.get_number_in_progress

```{autodoc2-docstring} eon.communicator.MPI.get_number_in_progress
```

````

````{py:method} cancel_state(state)
:canonical: eon.communicator.MPI.cancel_state

````

`````

`````{py:class} Local(scratchpath, client, ncpus, bundle_size, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.communicator.Local

Bases: {py:obj}`eon.communicator.Communicator`

```{autodoc2-docstring} eon.communicator.Local
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.communicator.Local.__init__
```

````{py:method} cleanup()
:canonical: eon.communicator.Local.cleanup

```{autodoc2-docstring} eon.communicator.Local.cleanup
```

````

````{py:method} get_results(resultspath, keep_result)
:canonical: eon.communicator.Local.get_results

```{autodoc2-docstring} eon.communicator.Local.get_results
```

````

````{py:method} check_job(job)
:canonical: eon.communicator.Local.check_job

```{autodoc2-docstring} eon.communicator.Local.check_job
```

````

````{py:method} submit_jobs(data, invariants)
:canonical: eon.communicator.Local.submit_jobs

```{autodoc2-docstring} eon.communicator.Local.submit_jobs
```

````

````{py:method} cancel_state(state)
:canonical: eon.communicator.Local.cancel_state

````

````{py:method} get_queue_size()
:canonical: eon.communicator.Local.get_queue_size

````

````{py:method} get_number_in_progress()
:canonical: eon.communicator.Local.get_number_in_progress

```{autodoc2-docstring} eon.communicator.Local.get_number_in_progress
```

````

`````

`````{py:class} Script(scratch_path, bundle_size, name_prefix, scripts_path, queued_jobs_cmd, cancel_job_cmd, submit_job_cmd, config: eon.config.ConfigClass = EON_CONFIG)
:canonical: eon.communicator.Script

Bases: {py:obj}`eon.communicator.Communicator`

```{autodoc2-docstring} eon.communicator.Script
```

```{rubric} Initialization
```

```{autodoc2-docstring} eon.communicator.Script.__init__
```

````{py:method} save_jobids()
:canonical: eon.communicator.Script.save_jobids

```{autodoc2-docstring} eon.communicator.Script.save_jobids
```

````

````{py:method} get_results(resultspath, keep_result)
:canonical: eon.communicator.Script.get_results

```{autodoc2-docstring} eon.communicator.Script.get_results
```

````

````{py:method} check_command(status, output, cmdname)
:canonical: eon.communicator.Script.check_command

```{autodoc2-docstring} eon.communicator.Script.check_command
```

````

````{py:method} submit_jobs(data, invariants)
:canonical: eon.communicator.Script.submit_jobs

````

````{py:method} cancel_state(state)
:canonical: eon.communicator.Script.cancel_state

````

````{py:method} get_queued_jobs()
:canonical: eon.communicator.Script.get_queued_jobs

```{autodoc2-docstring} eon.communicator.Script.get_queued_jobs
```

````

````{py:method} get_number_in_progress()
:canonical: eon.communicator.Script.get_number_in_progress

```{autodoc2-docstring} eon.communicator.Script.get_number_in_progress
```

````

````{py:method} get_queue_size()
:canonical: eon.communicator.Script.get_queue_size

````

`````
