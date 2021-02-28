#!/usr/bin/env nix-shell
#!nix-shell -i python3 -p "python38.withPackages(ps: with ps;[ numpy sh fire ])"

import textwrap
import re
import fire
import numpy as np
from pathlib import Path

def ams_stats(fname):
    cpu_re=re.compile('Total cpu time:\s*\d*.\d*')
    with open(Path(fname).parent/'ams_output') as f:
        amsc = f.read()
        allTimes = re.findall(cpu_re, amsc)
        allTimes = [float(x.split()[-1]) for x in allTimes]
    return allTimes

def print_stats(fname, ams=None):
    outstring = []
    with open(fname) as f:
        fcontents = f.read()
        totEcalls = fcontents.count('NORMAL TERMINATION')
        totDrops = fcontents.count('dropping')
        re_otime = re.compile('optimize time:\s*\d*.\d*')
        totTime = float(re.search(re.compile('real\s*\d*.\d*'), fcontents)[0].split()[-1])
        allOpt = []
        for line in re.findall(re_otime, fcontents):
            allOpt.append(float(line.split()[-1]))
        if (ams=="ams"):
            pSTime = np.array(ams_stats(fname))
            mSampleT = np.around(np.mean(pSTime)) # Rounded to ensure equality
            if mSampleT == 0:
                mSampleT = np.mean(pSTime) # No crashes
            overhead = totTime-(mSampleT*totEcalls)
            cost = (totTime-(mSampleT*totEcalls))/mSampleT
            outstring = f'''
            Average Single PES Sample: {np.round(mSampleT,3)}
            True PES Samples: {totEcalls}
            Dropped Samples: {totDrops}
            Hyperparameter Optimizations: {len(allOpt)}
            Hyperparameter Time Taken: {np.round(np.array(allOpt).sum(), 3)}
            Total Time (wall time): {totTime}
            Total PES Sample Time: {np.round(mSampleT*totEcalls, 3)}
            Wasted PES Sample Time: {np.round(mSampleT*totDrops, 3)}
            Overhead: {np.round(overhead, 3)}
            Cost: {np.round(cost,3)}
            '''
        else:
            outstring = f'''
            True PES Samples: {totEcalls}
            Dropped Samples: {totDrops}
            Hyperparameter Optimizations: {len(allOpt)}
            Hyperparameter Time Taken: {np.round(np.array(allOpt).sum(), 3)}
            Total Time (wall time): {totTime}
            '''
    return textwrap.dedent(outstring)

if __name__ == '__main__':
  fire.Fire(print_stats)
