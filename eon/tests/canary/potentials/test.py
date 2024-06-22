#!/usr/bin/env python
import shutil
import os
import sys

if not os.path.isfile("../../../client/client"):
    print("Client binary missing")
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
    shutil.copy(os.path.join('structs/',test['file']), 'pos.con')

    #write the config
    conf = open('config.ini','w')
    print('[Main]', file=conf)
    print('job = point', file=conf)
    print('[Potential]', file=conf)
    print('potential =', test['potential'], file=conf)
    conf.close()

    os.system("../../../client/client > /dev/null")

    result = open("results.dat",'r')
    energy = float(result.readline().split()[0])
    result.close()

    test['pass'] = abs(energy - test['energy']) <= test['tolerance']

allpassed = True
print("%-30s\t%-20s\t%-5s" % ("File", "Potential", "Pass"))
for test in tests:
    out = "%-30s\t%-20s\t%-5s" % (test['file'], test['potential'], test['pass'])
    if test['pass']:
        print('\033[92m'+out+'\033[0m')
    else:
        print('\033[91m'+out+'\033[0m')

    if not test['pass']:
        allpassed = False

os.unlink('config.ini')
os.unlink('pos.con')
os.unlink('results.dat')

if not allpassed:
    sys.exit(1)
