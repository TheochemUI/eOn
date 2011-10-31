#!/usr/bin/env python
import numpy

f = open('processtable')
processes = {}
for line in f:
    fields = line.split()
    if fields[0] == 'proc': continue
    number = int(fields[0])
    barrier = float(fields[6])
    frequency = int(fields[8])+1
    processes[number] = {'barrier':barrier, 'frequency':frequency} 

f.close()

f = open('search_results.txt')
f.readline()
f.readline()
for line in f:
    fields = line.split()
    result_string = fields[-1]
    if '-' not in result_string:
        continue
    first_part = result_string.split('-')[0]
    if first_part not in ("good", "repeat"):
        continue

    id = int(result_string.split('-')[1])
    saddle_forcecalls = int(fields[4])
    if 'fcs' not in processes[id]:
        processes[id]['fcs'] = []
    processes[id]['fcs'].append(saddle_forcecalls)

energies = []
for p in processes.itervalues():
    for i in range(p['frequency']):
        energies.append(p['barrier'])
    p['avg_fcs'] = float(sum(p['fcs'])) / p['frequency']
    p['stddev_fcs'] = numpy.std(p['fcs'])

bins = numpy.arange(min(energies)-0.05, max(energies)+.05, .05)
hist, bins = numpy.histogram(energies, bins=bins)
center = (bins[:-1]+bins[1:])/2.0

print 'set xlabel "Energy (eV)"'
print 'set ylabel "Frequncy"'
print 'set y2label "Force Calls"'
print 'set ytics nomirror'
print 'set y2tics'
print 'plot "-" with boxes title "Frequency'
for i in xrange(len(hist)):
    print "%.4f %i" % (center[i], hist[i])
