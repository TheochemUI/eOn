# Command Line Interface

```{versionchanged} 3.0
The CLI no longer handles jobs differently from the configuration file variant.
```

This documentation provides instructions for using the `eonclient` command line
interface (CLI).

## Basic Usage

```bash
eonclient [options]
```

## Options

- `-v, --version`
  **Description:** Print version information and exit.

- `-m, --minimize`
  **Description:** Run a minimization job using the `inputConfile` and save the
  results to `outputConfile`.

- `-s, --single`
  **Description:** Calculate the single point energy of the `inputConfile`.

- `-c, --compare`
  **Description:** Compare the structures between `inputConfile` and
  `outputConfile`.

- `--input=FILES`
  **Description:** Specify input files, separated by commas. Defaults to
  `pos.con` if not provided.

- `-o, --output=FILE`
  **Description:** Specify the output file.

- `--opt=METHOD`
  **Description:** Set the optimization method. Defaults to `cg`.

- `-f, --force=VALUE`
  **Description:** Set the convergence force for optimization. Default is
  `0.001`.

- `-t, --tolerance=VALUE`
  **Description:** Set the distance tolerance for structure comparison. Default
  is `0.1`.

- `--check_rotation`
  **Description:** Check rotation for structure comparison. Default is `true`.

- `-p, --potential=VALUE`
  **Description:** Specify the potential to use (e.g., `qsc`, `lj`, `eam_al`).

- `-q, --quiet`
  **Description:** Enable quiet mode. Suppresses version and timing information.

- `--conf=FILE`
  **Description:** Specify a configuration file to load. Default is
  `config.toml`.

- `-h, --help`
  **Description:** Print usage information and exit.

## Example Commands

### Running a Minimization Job

```{code-block} bash
eonclient --minimize --input=pos.con --output=output.con --potential=lj
```

### Running a Single Point Energy Calculation

```{code-block} bash
eonclient --single --input=pos.con --potential=qsc
```

### Comparing Structures

```{code-block} bash
eonclient --compare --input=pos.con,client/example/pos.con --output=output.con --check_rotation=false
```

### Running in Quiet Mode

```{code-block} bash
eonclient -s --input=pos.con --quiet --potential=lj
Energy:         -8.924553940551
Max atom force: 4.812074984293e+00
```

As compared to the timing information which is provided by default.

```{code-block} bash
eonclient -s --input=pos.con --potential=lj   
EON Client
VERSION: 1d84c297
BUILD DATE: Wed Aug 21 11:22:30 AM GMT 2024

OS: linux
Arch: x86_64
Hostname: rgx1gen11
PID: 1676485
DIR: /home/rgoswami/Git/Github/TheochemUI/eOn
Energy:         -8.924553940551
Max atom force: 4.812074984293e+00
timing information:
real      0.007 seconds
user      0.003 seconds
sys       0.010 seconds
```
