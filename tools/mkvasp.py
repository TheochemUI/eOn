#!/usr/bin/env python
import sys
import os

num_clients = int(sys.argv[1])
path = sys.argv[2]

for i in xrange(num_clients):
    dir_name = os.path.join(path, "vasp%3.3i"%i)
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
