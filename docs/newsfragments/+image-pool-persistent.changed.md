The NEB image pool keeps its threads for the life of the process and sizes itself to the cores the process may run on (the Slurm or taskset affinity mask) instead of every core of the node.
