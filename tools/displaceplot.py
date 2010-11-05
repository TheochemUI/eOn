#!/usr/bin/env python

import os

r = open("results.dat", 'r')
g = open("good.dat", 'w')
b = open("bad.dat", 'w')
lines = r.readlines()
data = []
for line in lines[1:]:
    split = line.split()
    if split[4] == 'good':
        g.write("%s    %s\n" % (split[0], split[1]))
    else:
        b.write("%s    %s\n" % (split[0], split[1]))
    
g.close()
b.close()
r.close()

os.system

dp = open('dp.plot', 'w')
dp.write('set xlabel "cutoff"\n')
dp.write('set ylabel "magnitude"\n')
dp.write('plot "bad.dat", "good.dat"\n')
dp.write('set term postscript eps color blacktext "Helvetica" 24\n')
dp.write('set output "displace.eps"\n')
dp.write('plot "bad.dat", "good.dat"\n')
dp.close()

os.system("gnuplot dp.plot -persist")




    
    
