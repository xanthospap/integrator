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
    '-s',
    '--saveas',
    metavar='SAVEAS',
    dest='save_as',
    default=None,
    required=False,
    help='Save the file using this file(name)')

parser.add_argument(
    '--scale-pos',
    metavar='SCALE_POSITION',
    dest='scale_pos',
    default=1e0,
    required=False,
    type=float,
    help='Scale plot')

parser.add_argument(
    '--scale-vel',
    metavar='SCALE_VELOCITY',
    dest='scale_vel',
    default=1e0,
    required=False,
    type=float,
    help='Scale plot')

parser.add_argument(
    '--no-diff',
    action='store_true',
    dest='nodif',
    required=False,
    help='Plot state results instead of differences.')

parser.add_argument(
    '--with-diff',
    action='store_true',
    dest='withdif',
    required=False,
    help='Plot state results and differences (same plot).')

#parser.add_argument(
#    '--plot-norm',
#    action='store_true',
#    dest='pnorm',
#    required=False,
#    help='Plot acceleration differences norm on another plot.')
#
#parser.add_argument(
#    '--quite',
#    action='store_true',
#    dest='quiet',
#    help='Do not show plot(s)')

if __name__ == "__main__":
    args = parser.parse_args()
    scale_pos = args.scale_pos
    scale_vel = args.scale_vel

    t = []; lx1=[]; ly1=[]; lz1=[];
    lx2=[]; ly2=[]; lz2=[];
    lvx1=[]; lvy1=[]; lvz1=[];
    lvx2=[]; lvy2=[]; lvz2=[];

    for line in sys.stdin:
        if line[0] != '#':
            try:
                sec, x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, vx2, vy2, vz2 = [ float(x) for x in line.split(' ') ]
                t.append(sec)
                lx1.append(x1); ly1.append(y1); lz1.append(z1);
                lx2.append(x2); ly2.append(y2); lz2.append(z2);
                lvx1.append(vx1); lvy1.append(vy1); lvz1.append(vz1);
                lvx2.append(vx2); lvy2.append(vy2); lvz2.append(vz2);
            except:
                pass
        else:
            if line.startswith('#title'):
                title = line.replace('#title','').strip()

    ## do nothing on empty input
    if len(t) <= 1: sys.exit(1)
    
    fig, axs = plt.subplots(3, 2, sharex=True)
    fig.subplots_adjust(hspace=0)
    group_labels = []

## Plot state results (not differences)
    if (args.nodif or args.withdif):
        print("Plotting state results ...")
        axs[0,0].plot(t, [scale_pos*z for z in lx1])
        axs[0,0].plot(t, [scale_pos*z for z in lx2])
        axs[0,1].plot(t, [scale_vel*z for z in lvx1])
        axs[0,1].plot(t, [scale_vel*z for z in lvx2])
        #
        axs[1,0].plot(t, [scale_pos*z for z in ly1])
        axs[1,0].plot(t, [scale_pos*z for z in ly2])
        axs[1,1].plot(t, [scale_vel*z for z in lvy1])
        axs[1,1].plot(t, [scale_vel*z for z in lvy2])
        #
        axs[2,0].plot(t, [scale_pos*z for z in lz1])
        axs[2,0].plot(t, [scale_pos*z for z in lz2])
        axs[2,1].plot(t, [scale_vel*z for z in lvz1])
        axs[2,1].plot(t, [scale_vel*z for z in lvz2])
        
        xtextpos = (t[0] + t[1]) / 2
        def ytextpos(scale, list): return scale * sum(z for z in list[0:5]) / 5
        
        axs[0,0].text(xtextpos,ytextpos(scale_pos, lx1),r'$X$ [m]')
        axs[1,0].text(xtextpos,ytextpos(scale_pos, ly1),r'$Y$ [m]')
        axs[2,0].text(xtextpos,ytextpos(scale_pos, lz1),r'$Z$ [m]')
        axs[0,1].text(xtextpos,ytextpos(scale_vel,lvx1),r'$V_{X}$ [m/sec]')
        axs[1,1].text(xtextpos,ytextpos(scale_vel,lvy1),r'$V_{Y}$ [m/sec]')
        axs[2,1].text(xtextpos,ytextpos(scale_vel,lvz1),r'$V_{Z}$ [m/sec]')

        # group_labels.append("acc. group 1")
        # group_labels.append("acc. group 2")

## Plot differences per component
    if (args.withdif or (not args.nodif)):
        print("Plotting state diffs ...")
        axs[0,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(lx1,lx2)])
        axs[1,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(ly1,ly2)])
        axs[2,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(lz1,lz2)])
        axs[0,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvx1,lvx2)])
        axs[1,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvy1,lvy2)])
        axs[2,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvz1,lvz2)])

        xtextpos = (t[0] + t[1]) / 2
        def ytextpos(scale, list1, list2):
            return scale*(sum(z for z in [ scale*(x[0]-x[1]) for x in zip(list1[0:5],list2[0:5]) ])) / 5
        
        axs[0,0].text(xtextpos,ytextpos(scale_pos, lx1, lx2),r'$X_{ref}-X$ [m]')
        axs[1,0].text(xtextpos,ytextpos(scale_pos, ly1, ly2),r'$Y_{ref}-Y$ [m]')
        axs[2,0].text(xtextpos,ytextpos(scale_pos, lz1, lz2),r'$Z_{ref}-Z$ [m]')
        axs[0,1].text(xtextpos,ytextpos(scale_vel,lvx1,lvx2),r'$V_{X_{ref}}-V_{X}$ [m/sec]')
        axs[1,1].text(xtextpos,ytextpos(scale_vel,lvy1,lvy2),r'$V_{Y_{ref}}-V_{Y}$ [m/sec]')
        axs[2,1].text(xtextpos,ytextpos(scale_vel,lvz1,lvz2),r'$V_{Z_{ref}}-V_{Z}$ [m/sec]')
        
        #plt.xlabel('Sec of Integration')
        #fig.text(0.0, 0.5, r'$\delta \ddot{r}_x$, $\delta \ddot{r}_y$ and $\delta \ddot{r}_z$ in $[m/sec^2]\times$'+'{:.1e}'.format(1/scale), va='center', rotation='vertical')
        
        # group_labels.append("acc. differences")

##  x-grids on
    #axs[0:].xaxis.grid(True, which='major')
    #axs[1:].xaxis.grid(True, which='major')
    #axs[2:].xaxis.grid(True, which='major')
    #axs[0:].xaxis.grid(True, which='minor')
    #axs[1:].xaxis.grid(True, which='minor')
    #axs[2:].xaxis.grid(True, which='minor')

    axs[2,0].set_xlabel("Sec of Integration (since t0) [sec]")
    axs[2,1].set_xlabel("Sec of Integration (since t0) [sec]")

    plt.show()

    if args.save_as:
        fn = os.path.join(args.save_as + 'raw_diffs')
        plt.savefig(fn+'.jpg')
