# Run a single eonclient job and capture results.dat as reference data

from pathlib import Path
import shutil

rule run_job:
    input:
        binary=EONCLIENT
    output:
        ref=str(OUTPUT_DIR / "{job}.dat")
    params:
        job_config=lambda wc: JOBS[wc.job],
        data_dir=DATA_DIR,
    run:
        job_cfg = params.job_config
        workdir = Path(f"workdirs/{wildcards.job}")
        workdir.mkdir(parents=True, exist_ok=True)

        data_source = job_cfg.get("data_source", "")
        source_dir = Path(str(params.data_dir)) / data_source

        # Copy all files from data source directory
        if source_dir.is_dir():
            for f in source_dir.iterdir():
                if f.is_file():
                    shutil.copy2(f, workdir / f.name)

        # Write config.ini if provided inline
        if "config" in job_cfg:
            (workdir / "config.ini").write_text(job_cfg["config"])

        # Rename pos file if specified (some jobs expect pos.con)
        pos_file = job_cfg.get("pos_file", "")
        if pos_file and (workdir / pos_file).exists():
            shutil.copy2(workdir / pos_file, workdir / "pos.con")

        # Run eonclient
        import subprocess
        result = subprocess.run(
            [str(Path(input.binary).resolve())],
            cwd=str(workdir),
            capture_output=True,
            text=True,
            timeout=120,
        )

        # Copy results.dat to output
        results_file = workdir / "results.dat"
        if results_file.exists():
            shutil.copy2(results_file, output.ref)
        else:
            # Some jobs (like FD) write to stdout instead
            Path(output.ref).write_text(result.stdout)
