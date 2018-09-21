#!/bin/sh
sbatch -A A-chgh -J "$1" -D "$2" -o ll_out -t "01:00:00" ~/code/eon/client/eonclient | awk '{print $4}'
