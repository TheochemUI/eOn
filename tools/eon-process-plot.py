#!/usr/bin/env python
import numpy
from sys import argv
from os.path import join

path = argv[1]

f = open(join(path, 'processtable'))
processes = {}
for line in f:
    fields = line.split()
    if fields[0] == 'proc': continue
    number = int(fields[0])
    barrier = float(fields[6])
    frequency = int(fields[8])+1
    processes[number] = {'barrier':barrier, 'frequency':frequency}

f.close()

f = open(join(path, 'search_results.txt'))
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
for p in processes.values():
    for i in range(p['frequency']):
        energies.append(p['barrier'])
    p['avg_fcs'] = float(sum(p['fcs'])) / p['frequency']
    p['stddev_fcs'] = numpy.std(p['fcs'])

bins = numpy.arange(min(energies)-0.05, max(energies)+.05, .05)
hist, bins = numpy.histogram(energies, bins=bins)
center = (bins[:-1]+bins[1:])/2.0

bin_map = {}
for id,p in processes.items():
    for i in range(len(bins)-1):
        if p['barrier'] >= bins[i] and p['barrier'] < bins[i+1]:
            if i not in bin_map:
                bin_map[i] = []
            bin_map[i].append(id)

print('set xlabel "Energy (eV)"')
print('set ylabel "Frequency"')
print('set y2label "Force Calls"')
print('set ytics nomirror')
print('set y2tics')
print('set autoscale y')
print('set autoscale y2')
print('set boxwidth 0.05')
print('plot "-" with boxes fs solid 0.5 lt 3 title "Frequency" axes x1y1, "-" with yerrorbars title "Force Calls" axes x1y2 lt -1')
for i in range(len(hist)):
    if i not in bin_map:
        continue

    print("%.4f %i" % (center[i], hist[i]))

print("e")

for i in range(len(hist)):
    if i not in bin_map:
        continue
    fcs = []
    for id in bin_map[i]:
       fcs.extend(processes[id]['fcs'])
    print("%.4f %.4f %.4f" % (center[i], numpy.average(fcs), numpy.std(fcs)))
