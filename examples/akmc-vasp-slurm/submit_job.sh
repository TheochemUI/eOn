#!/bin/bash
sbatch \
	--job-name="$1" \
	--output=eonclient_stdout.log \
	--qos=debug \
	--nodes=1 \
	--time=0:10:00 \
	-C haswell \
	--account=account_name \
	--mail-type=ALL \
	--mail-user=your@email \
	--chdir="$2" \
	client_script.sh | awk '{print $4}'
#	~/software/eon/bin/eonclient | awk '{print $4}'
#	module load gcc; module load vasp/5.4.1_vtst-gcc; module load ase; ~/software/eon/bin/eonclient | awk '{print $4}'
#	module load python/2.7-anaconda-5.2 vasp/5.4.1_vtst ase
