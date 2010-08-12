#!/usr/bin/env awk -f

# generates a gnuplot script and data for energy versus time
# use by running from the simulation directory as such:
# evt.awk states/state_table dynamics.txt | gnuplot

BEGIN { 
    print "set title \"Energy of Minima versus Time\""
    print "set xlabel \"time(s)\""
    print "set ylabel \"Energy (eV)\""
    print "set nokey"
    print "plot \"-\" with lines" 
}

$3 !~ /[0-9]+/ {e[$1]=$2} 
$3 ~ /[0-9]+/ {print $3, e[$1];}
