#!/bin/bash
set -e

cp config_passed_dimer_improved_lbfgs.ini config_passed.ini
../../client/client > dimer_lbfgs.out
cp mode.dat mode_dimer_improved_lbfgs.dat

cp config_passed_dimer_improved_cg.ini config_passed.ini
../../client/client > dimer_cg.out
cp mode.dat mode_dimer_improved_cg.dat

cp config_passed_dimer_improved_sd.ini config_passed.ini
../../client/client > dimer_sd.out
cp mode.dat mode_dimer_improved_sd.dat

cp config_passed_lanczos.ini config_passed.ini
../../client/client > lanczos.out
cp mode.dat mode_lanczos.dat

sd_lanczos=`../../tools/modedot.py mode_dimer_improved_sd.dat mode_lanczos.dat`
cg_lanczos=`../../tools/modedot.py mode_dimer_improved_cg.dat mode_lanczos.dat`
lbfgs_lanczos=`../../tools/modedot.py mode_dimer_improved_lbfgs.dat mode_lanczos.dat`

echo "Dot product of converged mode:"
printf "dimer sd    vs lanczos: %.3f\n" $sd_lanczos
printf "dimer cg    vs lanczos: %.3f\n" $cg_lanczos
printf "dimer lbfgs vs lanczos: %.3f\n" $lbfgs_lanczos
