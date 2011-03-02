#!/usr/bin/env python
import shutil
import os
import subprocess
import sys

if not os.path.isfile("../../client/client"):
    print "Client binary missing"
    sys.exit(1)
tab = open('energies')

tests = []
for line in tab:
    line = line.split()
    # 0 - filename
    # 1 - potential
    # 2 - tolerance
    tests.append({'file':line[0], 'potential':line[1], 'energy':float(line[2]), 
        'tolerance':float(line[3])})
tab.close()

for test in tests:
    #copy the reactant
    shutil.copy(os.path.join('structs/',test['file']), 'reactant_passed.con')    
    
    #write the config
    conf = open('config_passed.ini','w')
    print >> conf, '[Main]'
    print >> conf, 'job = point'
    print >> conf, 'potential =', test['potential']
    conf.close()
    
    os.system("../../client/client > /dev/null")

    result = open("results.dat",'r')
    energy = float(result.readline().split()[1])
    result.close()

    test['pass'] = abs(energy - test['energy']) <= test['tolerance']

allpassed = True
print "%-30s\t%-20s\t%-5s" % ("File", "Potential", "Pass")
for test in tests:
    print "%-30s\t%-20s\t%-5s" % (test['file'], test['potential'], test['pass'])
    if not test['pass']:
        allpassed = False

os.unlink('config_passed.ini')
os.unlink('reactant_passed.con')
os.unlink('results.dat')

if not allpassed:
    sys.exit(1)

