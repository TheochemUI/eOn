In-process jobs return a ``JobResult`` (scalars plus ConFrame saddle and
product). ``results.dat`` is only the ``job_result_to_results_dat`` adapter
for cluster and HPC workers. ``parse_results`` still reads legacy files.
