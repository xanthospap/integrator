#! /usr/bin/python

import os
import sys
import math
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import argparse

parser = argparse.ArgumentParser(
    description=
    'Plot state integration diffs w.r.t reference results.',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Apr, 2024'''))

parser.add_argument(
    '-c',
    '--costg-ref',
    metavar='COSTG_REFERENCE_FILE',
    dest='costg_ref',
    default=None,
    required=True,
    help='The COSTG-benchmark reference file to check against.')

def parse_costg_fn(fn):
    lmjd=[]; lax1=[]; lay1=[]; laz1=[];
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            try:
                #print(len(line.split()))
                mjd, ax, ay, az = [ float(x) for x in line.split() ]
                #print("Done[1]");
                lmjd.append(mjd)
                lax1.append(ax); lay1.append(ay); laz1.append(az);
            except:
                print("Skipped line {:}".format(line.strip()))
                #pass
    return lmjd, lax1, lay1, laz1

if __name__ == "__main__":
    args = parser.parse_args()

    lmjd1=[]; lax1=[]; lay1=[]; laz1=[];
    for line in sys.stdin:
        if line[0] != '#':
            try:
                mjd, ax, ay, az = [ float(x) for x in line.split(' ') ]
                lmjd1.append(mjd)
                lax1.append(ax); lay1.append(ay); laz1.append(az);
            except:
                pass
        else:
            if line.startswith('#title'):
                title = line.replace('#title','').strip()


    fig, axs = plt.subplots(4, 1, sharex=True)
    fig.subplots_adjust(hspace=0)

    lmjd2, lax2, lay2, laz2 = parse_costg_fn(args.costg_ref)
    print("Read {:} records of from ref file".format(len(lmjd2)))
    assert(len(laz2) >= len(laz1))
    sz = len(lmjd1)
    axs[0].plot(lmjd1, [(z[0]-z[1]) for z in zip(lax2[0:sz],lax1)])
    axs[1].plot(lmjd1, [(z[0]-z[1]) for z in zip(lay2[0:sz],lay1)])
    axs[2].plot(lmjd1, [(z[0]-z[1]) for z in zip(laz2[0:sz],laz1)])
    axs[3].plot(lmjd1, [(z[0]-z[1])*86400e0 for z in zip(lmjd2[0:sz],lmjd1)])
    
    plt.show()
