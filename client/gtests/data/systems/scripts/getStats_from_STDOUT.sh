#!/usr/bin/env nix-shell
#!nix-shell -i python3 -p "python38.withPackages(ps: with ps;[ numpy sh fire ])"

import textwrap
import re
import fire
import numpy as np
from pathlib import Path

def print_stats(fname, ams=None):
    outstring = []
    with open(fname) as f:
        fcontents = f.read()
        totEcalls = fcontents.count('NORMAL TERMINATION')
        totDrops = fcontents.count('dropping')
        re_otime = re.compile('optimize time: \d*.\d*')
        totTime = float(re.search(re.compile('real    \d*.\d*'), fcontents)[0].split()[-1])
        allOpt = []
        for line in re.findall(re_otime, fcontents):
            allOpt.append(float(line.split()[-1]))
        if (ams=="ams"):
            with open(Path(fname).parent/'ams_output') as acontent:
                aconts = acontent.read()
                pSTime = float(re.search(re.compile('Total cpu time:\s*\d*.\d*'), aconts)[0].split()[-1])
                outstring = f'''
                Single PES Sample: {pSTime}
                True PES Samples: {totEcalls}
                Dropped Samples: {totDrops}
                Hyperparameter Optimizations: {len(allOpt)}
                Hyperparameter Time Taken: {np.array(allOpt).sum()}
                Total Time (wall time): {totTime}
                Total PES Sample Time: {pSTime*totEcalls}
                '''
        else:
            outstring = f'''
            True PES Samples: {totEcalls}
            Dropped Samples: {totDrops}
            Hyperparameter Optimizations: {len(allOpt)}
            Hyperparameter Time Taken: {np.array(allOpt).sum()}
            Total Time (wall time): {totTime}
            '''
    return textwrap.dedent(outstring)

if __name__ == '__main__':
  fire.Fire(print_stats)
