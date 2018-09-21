sbatch -A A-chgh -J "$1" -D "$2" -n 48 -N 1 -p skx-dev -o ll_out -t "01:00:00" --wrap="~/code/eon/client/eonclient" > sbatch.out
cat sbatch.out | tail -1 | awk '{print $4}'
