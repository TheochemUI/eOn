#!/usr/local/Cellar/gnuplot/5.4.2/bin/gnuplot --persist

set terminal postscript enhanced color 
set output "neb_paths.eps"
set style line 1 linewidth 2.000 dashtype solid pointtype 7 pointsize 1.500 pointinterval -1
set encoding utf8
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit autofreq  font ",10"
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq  font ",10"
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
set title "NEB-GPR path evolution" 
set title  font "" textcolor lt -1 norotate
set xlabel "path (Å)" 
set xlabel  font ",16" textcolor lt -1 norotate
set ylabel "Energy (eV)" 
set ylabel  font ",16" textcolor lt -1 rotate
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
#unset key
GNUTERM = "qt"

filename(n) = sprintf("neb_gpr_%d", n)
plot for [i=1:ARG1] filename(i).'.dat' using 2:3 w lp ls 1 lc i title 'Path: '.i
