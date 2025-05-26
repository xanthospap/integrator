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
    default=1e3,
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

parser.add_argument(
    '--compare-to',
    metavar='COMPARE_TO_FILE',
    dest='compare_to',
    required=False,
    default=None,
    help='Plot results versus a previous run.')

parser.add_argument(
    '-f',
    '--from-file',
    metavar='READ_FROM_FILE',
    dest='from_file',
    required=False,
    default=None,
    help='Read data to plot from file instead of STDIN.')

parser.add_argument(
    '--xhours',
    action='store_true',
    dest='xaxh',
    help='Make x-axis hours instead of seconds.')

defaultPlotOptions = {
    "style_sheet": "dark_background",
    "data_points_color": "blue",
    "data_points_color2": "black",
    "data_points_line_width": 1.2,
    "line_color": "blue",
    "line_color2": "black",
    "line_style": "--",
    # "error_bar_color": "#ff9999",
    "error_bar_color": (1.0, 0.3, 0.3, 0.3),
    "error_bar_color2": (1.0, 0.3, 0.6, 0.6),
    "error_bar_width": 2,
    "error_bar_capsize": 3,
    "shaded_error_bars": True,
    "error_bar_alpha": 0.2,
}

def parse(source, sec2hr=False):
    t = []; lx1=[]; ly1=[]; lz1=[];
    lx2=[]; ly2=[]; lz2=[];
    lvx1=[]; lvy1=[]; lvz1=[];
    lvx2=[]; lvy2=[]; lvz2=[];
    mean_height=0e0; mean_velo=0e0;
    N = 0

    for line in source:
        if line[0] != '#':
            try:
                sec, x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, vx2, vy2, vz2 = [ float(x) for x in line.split(' ') ]
                t.append(sec if sec2hr is False else sec/3600.)
                lx1.append(x1); ly1.append(y1); lz1.append(z1);
                lx2.append(x2); ly2.append(y2); lz2.append(z2);
                lvx1.append(vx1); lvy1.append(vy1); lvz1.append(vz1);
                lvx2.append(vx2); lvy2.append(vy2); lvz2.append(vz2);
                N += 1
                mean_height = mean_height * (N-1)/N + np.linalg.norm(np.array((x1,y1,z1))*1e-3)/N
                mean_velo = mean_velo * (N-1)/N + np.linalg.norm(np.array((vx1,vy1,vz1)))/N
            except:
                pass
        else:
            if line.startswith('#title'):
                title = line.replace('#title','').strip()
    return t, lx1, ly1, lz1, lvx1, lvy1, lvz1, lx2, ly2, lz2, lvx2, lvy2, lvz2, mean_height, mean_velo

if __name__ == "__main__":
    args = parser.parse_args()
    scale_pos = args.scale_pos
    scale_vel = args.scale_vel

    plotOptions = defaultPlotOptions

    # read from file or stdin
    if args.from_file is not None:
        with open(args.from_file, 'r') as f:
            t, lx1, ly1, lz1, lvx1, lvy1, lvz1, lx2, ly2, lz2, lvx2, lvy2, lvz2, mean_height, mean_velo = parse(f, args.xaxh)
    else:
        t, lx1, ly1, lz1, lvx1, lvy1, lvz1, lx2, ly2, lz2, lvx2, lvy2, lvz2, mean_height, mean_velo = parse(sys.stdin, args.xaxh)

    ## do nothing on empty input
    if len(t) <= 1: sys.exit(1)
    
    fig, axs = plt.subplots(3, 2, sharex=True)
    fig.subplots_adjust(hspace=0)
    group_labels = []

## Plot results compared to another file
    if args.compare_to is not None:
        st, slx1, sly1, slz1, slvx1, slvy1, slvz1, slx2, sly2, slz2, slvx2, slvy2, slvz2, _, _ = parse(args.compare_to, args.xaxh)
        axs[0,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(lx1,lx2)]  ,label=r'$\delta x$ current')
        axs[1,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(ly1,ly2)]  ,label=r'$\delta y$ current')
        axs[2,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(lz1,lz2)]  ,label=r'$\delta z$ current')
        axs[0,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvx1,lvx2)],label=r'$\delta v_x$ current')
        axs[1,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvy1,lvy2)],label=r'$\delta v_y$ current')
        axs[2,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvz1,lvz2)],label=r'$\delta v_z$ current')
        
        axs[0,0].plot(st, [scale_pos*(z[0]-z[1]) for z in zip(slx1, slx2)] ,label=r'$\delta x$ {:}'.format(args.compare_to))
        axs[1,0].plot(st, [scale_pos*(z[0]-z[1]) for z in zip(sly1, sly2)] ,label=r'$\delta y$ {:}'.format(args.compare_to))
        axs[2,0].plot(st, [scale_pos*(z[0]-z[1]) for z in zip(slz1, slz2)] ,label=r'$\delta z$ {:}'.format(args.compare_to))
        axs[0,1].plot(st, [scale_vel*(z[0]-z[1]) for z in zip(slvx1,slvx2)],label=r'$\delta v_x$ {:}'.format(args.compare_to))
        axs[1,1].plot(st, [scale_vel*(z[0]-z[1]) for z in zip(slvy1,slvy2)],label=r'$\delta v_y$ {:}'.format(args.compare_to))
        axs[2,1].plot(st, [scale_vel*(z[0]-z[1]) for z in zip(slvz1,slvz2)],label=r'$\delta v_z$ {:}'.format(args.compare_to))
        
        axs[0,0].legend(loc='upper left')
        axs[1,0].legend(loc='upper left')
        axs[2,0].legend(loc='upper left')
        axs[0,1].legend(loc='upper left')
        axs[1,1].legend(loc='upper left')
        axs[2,1].legend(loc='upper left')

## Plot state results (not differences)
    if (args.nodif or args.withdif) and not args.compare_to:
        print("Plotting state results ...")
        axs[0,0].plot(t, [scale_pos*z for z in lx1] ,color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[0,0].plot(t, [scale_pos*z for z in lx2] ,color=plotOptions["line_color2"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[0,1].plot(t, [scale_vel*z for z in lvx1],color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[0,1].plot(t, [scale_vel*z for z in lvx2],color=plotOptions["line_color2"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[1,0].plot(t, [scale_pos*z for z in ly1] ,color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[1,0].plot(t, [scale_pos*z for z in ly2] ,color=plotOptions["line_color2"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[1,1].plot(t, [scale_vel*z for z in lvy1],color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[1,1].plot(t, [scale_vel*z for z in lvy2],color=plotOptions["line_color2"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[2,0].plot(t, [scale_pos*z for z in lz1] ,color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1,label="group A")
        axs[2,0].plot(t, [scale_pos*z for z in lz2] ,color=plotOptions["line_color2"],linestyle=plotOptions["line_style"],zorder=1,label="group B")
        axs[2,1].plot(t, [scale_vel*z for z in lvz1],color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        axs[2,1].plot(t, [scale_vel*z for z in lvz2],color=plotOptions["line_color2"],linestyle=plotOptions["line_style"],zorder=1,label="_nolegend_")
        
        axs[0,0].scatter(t, [scale_pos*z for z in lx1] ,facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[0,0].scatter(t, [scale_pos*z for z in lx2] ,facecolors=plotOptions['data_points_color2'],edgecolors=plotOptions['error_bar_color2'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[0,1].scatter(t, [scale_vel*z for z in lvx1],facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[0,1].scatter(t, [scale_vel*z for z in lvx2],facecolors=plotOptions['data_points_color2'],edgecolors=plotOptions['error_bar_color2'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[1,0].scatter(t, [scale_pos*z for z in ly1] ,facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[1,0].scatter(t, [scale_pos*z for z in ly2] ,facecolors=plotOptions['data_points_color2'],edgecolors=plotOptions['error_bar_color2'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[1,1].scatter(t, [scale_vel*z for z in lvy1],facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[1,1].scatter(t, [scale_vel*z for z in lvy2],facecolors=plotOptions['data_points_color2'],edgecolors=plotOptions['error_bar_color2'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[2,0].scatter(t, [scale_pos*z for z in lz1] ,facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[2,0].scatter(t, [scale_pos*z for z in lz2] ,facecolors=plotOptions['data_points_color2'],edgecolors=plotOptions['error_bar_color2'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[2,1].scatter(t, [scale_vel*z for z in lvz1],facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        axs[2,1].scatter(t, [scale_vel*z for z in lvz2],facecolors=plotOptions['data_points_color2'],edgecolors=plotOptions['error_bar_color2'],linewidth=plotOptions['data_points_line_width'],zorder=2,label="_nolegend_")
        
        axs[2,0].legend(loc='upper left')
        axs[0,0].text(.01,.99,r'$X$'        ,transform=axs[0,0].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[1,0].text(.01,.99,r'$Y$'        ,transform=axs[1,0].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[2,0].text(.01,.99,r'$Z$'        ,transform=axs[2,0].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[0,1].text(.01,.99,r'$V_{X}$',transform=axs[0,1].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[1,1].text(.01,.99,r'$V_{Y}$',transform=axs[1,1].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[2,1].text(.01,.99,r'$V_{Z}$',transform=axs[2,1].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)

        # Add invisible axes for left and right column labels
        fig.text(0.04, 0.5, '[m]', va='center', rotation='vertical')
        fig.text(0.96, 0.5, '[mm/sec]', va='center', rotation='vertical')
        
## Plot differences per component
    if (args.withdif or (not args.nodif)) and not args.compare_to:
        print("Plotting state diffs ...")
        axs[0,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(lx1,lx2)]  ,color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1)
        axs[1,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(ly1,ly2)]  ,color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1)
        axs[2,0].plot(t, [scale_pos*(z[0]-z[1]) for z in zip(lz1,lz2)]  ,color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1)
        axs[0,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvx1,lvx2)],color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1)
        axs[1,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvy1,lvy2)],color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1)
        axs[2,1].plot(t, [scale_vel*(z[0]-z[1]) for z in zip(lvz1,lvz2)],color=plotOptions["line_color"],linestyle=plotOptions["line_style"],zorder=1)
        
        axs[0,0].scatter(t, [scale_pos*(z[0]-z[1]) for z in zip(lx1,lx2)]  ,facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2)
        axs[1,0].scatter(t, [scale_pos*(z[0]-z[1]) for z in zip(ly1,ly2)]  ,facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2)
        axs[2,0].scatter(t, [scale_pos*(z[0]-z[1]) for z in zip(lz1,lz2)]  ,facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2)
        axs[0,1].scatter(t, [scale_vel*(z[0]-z[1]) for z in zip(lvx1,lvx2)],facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2)
        axs[1,1].scatter(t, [scale_vel*(z[0]-z[1]) for z in zip(lvy1,lvy2)],facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2)
        axs[2,1].scatter(t, [scale_vel*(z[0]-z[1]) for z in zip(lvz1,lvz2)],facecolors=plotOptions['data_points_color'],edgecolors=plotOptions['error_bar_color'],linewidth=plotOptions['data_points_line_width'],zorder=2)

        axs[0,0].text(.01,.99,r'$X_{ref}-X$'        ,transform=axs[0,0].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[1,0].text(.01,.99,r'$Y_{ref}-Y$'        ,transform=axs[1,0].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[2,0].text(.01,.99,r'$Z_{ref}-Z$'        ,transform=axs[2,0].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[0,1].text(.01,.99,r'$V_{X_{ref}}-V_{X}$',transform=axs[0,1].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[1,1].text(.01,.99,r'$V_{Y_{ref}}-V_{Y}$',transform=axs[1,1].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)
        axs[2,1].text(.01,.99,r'$V_{Z_{ref}}-V_{Z}$',transform=axs[2,1].transAxes,verticalalignment='top', horizontalalignment='left',fontsize=10,clip_on=True)

        # Add invisible axes for left and right column labels
        fig.text(0.04, 0.5, '[m]', va='center', rotation='vertical')
        fig.text(0.96, 0.5, '[mm/sec]', va='center', rotation='vertical')

        for ax in axs[:, 1]:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position('right')

    for ax in axs.flat: ax.grid(True)
    tu = 'Seconds' if args.xaxh is False else 'Hours'
    axs[2,0].set_xlabel(f"{tu} of Integration (since t0)")
    axs[2,1].set_xlabel(f"{tu} of Integration (since t0)")
    fig.suptitle('Mean Altitude is {:.1f}km, mean velocity is {:.1f}m/sec'.format(mean_height-6378.1, mean_velo), fontsize=16)
    # Left column (column 0)
    left_axes = axs[:, 0]
# Right column (column 1)
    right_axes = axs[:, 1]
# Compute y-limits for each column
    left_ymins, left_ymaxs = zip(*[ax.get_ylim() for ax in left_axes])
    right_ymins, right_ymaxs = zip(*[ax.get_ylim() for ax in right_axes])
# Get global min and max for each column
    left_ylim = (min(left_ymins), max(left_ymaxs))
    right_ylim = (min(right_ymins), max(right_ymaxs))
# Apply y-limits
    for ax in left_axes:
        ax.set_ylim(left_ylim)
    for ax in right_axes:
        ax.set_ylim(right_ylim)

    # reduce width spacing between columns
    fig.subplots_adjust(wspace=0.1)

    plt.show()

    if args.save_as:
        fn = os.path.join(args.save_as + 'raw_diffs')
        plt.savefig(fn+'.jpg')
