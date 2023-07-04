#!/usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
from optparse import OptionParser
import configparser
import subprocess

def initialize(rads, mags, searches):
    os.mkdir('original')
    things = glob.glob('*')
    for thing in things:
        if thing != 'original' and not thing.startswith('samples-'):
            shutil.move(thing, 'original')
    os.mkdir('samples')
    write_metadata(0, rads, 0, mags, searches)

def restore():
    if not os.path.exists('original'):
        print('Nothing to restore!')
        return
    os.mkdir('samples/continue')
    things = glob.glob('*')
    for thing in things:
        if thing not in ['original', 'samples'] and not thing.startswith('samples-'):
            shutil.move(thing, 'samples/continue')
    things = glob.glob('original/*')
    for thing in things:
        shutil.move(thing, '.')
    os.rmdir('original')
    shutil.move('samples', 'samples-' + str(datetime.datetime.now()).replace(' ', '-'))

def read_metadata():
    f = open('sample-status.txt', 'r')
    lines = f.readlines()
    f.close()
    radius_index = int(lines[0].split()[2])
    radii = [float(r) for r in lines[0].split()[4:]]
    magnitude_index = int(lines[1].split()[2])
    magnitudes = [float(m) for m in lines[1].split()[4:]]
    searches = int(lines[2].split()[1])
    return radius_index, radii, magnitude_index, magnitudes, searches

def write_metadata(radi, rads, magi, mags, searches):
    ss = open('sample-status.txt', 'w')
    ss.write('radii: index %d of %s\n' %(radi, ' '.join([str(r) for r in rads])))
    ss.write('magnitudes: index %d of %s\n' % (magi, ' '.join([str(m) for m in mags])))
    ss.write('searches: %d\n' % searches)
    ss.close()

def inc_params():
    radi, rads, magi, mags, s = read_metadata()
    radi += 1
    if radi >= len(rads):
        radi = 0
        magi += 1
    write_metadata(radi, rads, magi, mags, s)

def store_sample(location):
    os.mkdir(location)
    things = glob.glob('*')
    for thing in things:
        if thing not in ['original', 'samples', 'sample-status.txt'] and not thing.startswith('samples-'):
            shutil.move(thing, location)

def print_table():
    samples = glob.glob('samples/*')
    print('%10s  %10s  %10s  %s' % ('radius', 'magnitude', '#procs', '#goodsearch'))
    max_procs = 0
    for sample in samples:
        if not sample.startswith('samples/sample-'):
            continue
        max_procs = max(len(open('%s/states/0/processtable' % sample, 'r').readlines()) - 1, max_procs)
    for sample in samples:
        if not sample.startswith('samples/sample-'):
            continue
        rad = float(sample.split('-')[1])
        mag = float(sample.split('-')[2])
        total_searches = len(open('%s/states/0/search_results.txt' % sample, 'r').readlines()) - 2
        good_searches = sum([int(l.split()[8])+1 for l in open('%s/states/0/processtable' % sample, 'r').readlines()[1:]])
        procs = len(open('%s/states/0/processtable' % sample, 'r').readlines()) - 1
        print('%10.5f  %10.5f  %10d  %d' % (rad, mag, procs, good_searches))



def main():

    parser = OptionParser()
    parser.add_option('--restore', dest='restore', action='store_true',
                      help='restore your simulation files and cache the current sampling results',
                      default=False)
    parser.add_option('--radii', dest='radii', default='0.0,4.0,8.0',
                      action='store', type='string',
                      help='the displacement radii to sample, separated with commas and no spaces, e.g., --radius=0.0,4.0,8.0')
    parser.add_option('--magnitudes', dest='magnitudes', default='0.1,0.2,0.4',
                      action='store', type='string',
                      help='the displacement magnitudes to sample, separated with commas and no spaces, e.g., --magnitude=0.1,0.2,0.4')
    parser.add_option('--searches', dest='searches', default=50,
                      action='store', type=int,
                      help='the number of searches to perform per sample')
    (options, args) = parser.parse_args()

    if options.restore:
        restore()
        return

    if not os.path.exists('original'):
        if not os.path.exists('config.ini') or not os.path.exists('reactant.con'):
            print('This does not appear to be an eon project directory. Either config.ini or reactant.con or both were not found. Aborting.')
            sys.exit()
        rads = [float(r) for r in options.radii.split(',')]
        mags = [float(m) for m in options.magnitudes.split(',')]
        searches = int(options.searches)
        initialize(rads, mags, searches)
        print('Initialized new akmc displacement sampling with radii %s, magnitudes %s, and performing %s searches per sample.' % (options.radii, options.magnitudes, options.searches))
        return

    radi, rads, magi, mags, searches = read_metadata()

    # Make sure we're not finished sampling before continuing.
    if magi >= len(mags):
        print_table()
        print('Sampling complete.')
        sys.exit(1)
        return

    rad = rads[radi]
    mag = mags[magi]

    if not os.path.exists('config.ini') or not os.path.exists('reactant.con'):
        # Initialize the next sample.
        for thing in glob.glob('original/*'):
            shutil.copy(thing, '.')
        print('Sampling with radius %f and magnitude %f' % (rad, mag))
        parser = configparser.SafeConfigParser()
        parser.read('config.ini')
        parser.set('Saddle Search', 'displace_radius', str(rad))
        parser.set('Saddle Search', 'displace_magnitude', str(mag))
        parser.set('AKMC', 'confidence', '1.1')
        fob = open('config.ini', 'w')
        parser.write(fob)
        fob.close()

    # Register any available results.
    throwaway = subprocess.getoutput('eon -n')

    # Check the number of completed searches.
    n_completed_searches = 0
    if os.path.exists('states/0/search_results.txt'):
        n_completed_searches = len(open('states/0/search_results.txt').readlines()) - 2
    print('%d searches have been completed for this sample.' % n_completed_searches)
    if n_completed_searches >= searches:
        print('Sample finished.')
        store_sample('samples/sample-%f-%f' % (rad, mag))
        # Increment the sample parameters
        inc_attempt = inc_params()
        return

    print('Executing eon...')
    throwaway = subprocess.getoutput('eon')
    print('Done.')

if __name__ == '__main__':
    main()
