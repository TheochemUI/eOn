# client executable must be copied to same directory
# created data is stored in the subdirectory performance

generateData = 1
sampleSize = 500

import os
import shutil


f_positive = []
f_negative = []

def getLinesInFile(fileName):
    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()
    return lines

if generateData == 1:
    os.mkdir("data")
    for i in range(sampleSize):
        runClient = "./client"
        os.system(runClient)
        os.mkdir("data/"+str(i))
        shutil.move("results.dat","data/"+str(i))

for i in range(sampleSize):
    data = getLinesInFile("./data/"+str(i)+"/results.dat")
    f_positive.append(int(data[5].rsplit()[0]))
    f_negative.append(int(data[6].rsplit()[0]))

print "Average nr of force call per saddle point: " + str((sum(f_negative) + sum(f_positive)) /float(sampleSize))
