#! /usr/bin/python

import os
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import argparse

parser = argparse.ArgumentParser(
    description="Plot state integration diffs w.r.t reference results.",
    epilog=(
        """National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Apr, 2024"""
    ),
)

parser.add_argument(
    "-d",
    "--costg-data-dir",
    metavar="COSTG_DATA_DIR",
    dest="costg_dir",
    default=None,
    required=True,
    help="The COSTG-benchmark reference data folder to check against.",
)


def parse_costg_fn(fn):
    lmjd = []
    lax1 = []
    lay1 = []
    laz1 = []
    with open(fn, "r") as fin:
        for line in fin.readlines():
            try:
                mjd, ax, ay, az = [float(x) for x in line.split()]
                lmjd.append(mjd)
                lax1.append(ax)
                lay1.append(ay)
                laz1.append(az)
            except:
                # print("Skipped line {:}".format(line.strip()))
                pass
    return lmjd, lax1, lay1, laz1


def parse_acc():
    t = []
    a1 = []
    a2 = []
    a3 = []
    a4 = []
    a5 = []
    a6 = []
    a7 = []
    a8 = []
    a9 = []
    a10 = []
    for line in sys.stdin:
        if line[0] != "#":
            l = line.split()
            assert len(l) >= 25
            t.append(float(l[0]))
            a1.append(np.array([float(x) for x in l[1:4]]))
            a2.append(np.array([float(x) for x in l[4:7]]))
            a3.append(np.array([float(x) for x in l[7:10]]))
            a4.append(np.array([float(x) for x in l[10:13]]))
            a5.append(np.array([float(x) for x in l[13:16]]))
            a6.append(np.array([float(x) for x in l[16:19]]))
            a7.append(np.array([float(x) for x in l[19:22]]))
            a8.append(np.array([float(x) for x in l[22:25]]))
            a9.append(np.array([float(x) for x in l[25:28]]))
            a10.append(np.array([float(x) for x in l[28:31]]))
    return t, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10


def findCloseMatches(t1, t2, tol=1e-12):
    matches = []
    i, j = 0, 0
    n1, n2 = len(t1), len(t2)

    while i < n1 and j < n2:
        diff = t1[i] - t2[j]
        if abs(diff) < tol:
            matches.append((i, j))
            i += 1
            j += 1
        elif diff < 0:
            i += 1
        else:
            j += 1

    return matches


def genericData(costg_fn, ttest, ytest):
    lmjd, ax, ay, az = parse_costg_fn(costg_fn)
    t = []
    yr = []
    yt = []
    m = findCloseMatches(lmjd, ttest)
    for tpl in m:
        i = tpl[0]
        j = tpl[1]
        t.append(lmjd[i])
        yr.append(np.array([ax[i], ay[i], az[i]]))
        yt.append(ytest[j])
    print(f"Matched epochs: {len(t)}")
    return t, yr, yt


def genericPlotter(
    costg_fn, ttest, ytest, y_axis_label, title, plot_diff=True, optionsFn=None
):
    dps = 0.5
    mt, mref, mtst = genericData(costg_fn, ttest, ytest)
    fig, ax = plt.subplots(3, 1, sharex=True)
    if plot_diff:
        for i, c in zip(range(3), ["X", "Y", "Z"]):
            ax[i].scatter(mt, [a[i] - b[i] for a, b in zip(mtst, mref)], s=dps)
            ax[i].grid(True)
            ax[i].set_ylabel(f"{c}: {y_axis_label}")
    else:
        for i, c in zip(range(3), ["X", "Y", "Z"]):
            ax[i].scatter(mt, [a[i] for a in mref], s=dps, label="costg")
            ax[i].scatter(mt, [a[i] for a in mtst], s=dps, label="test")
            ax[i].grid(True)
            ax[i].set_ylabel(f"{c}: {y_axis_label}")
        ax[2].legend(loc="upper right")
    ax[2].set_xlabel("Epoch")
    ax[0].set_title(f"{title}")
    fig.subplots_adjust(hspace=0)
    return fig


costg_data_files = {
    "earth_gravity": {"ref": "02gravityfield_icrf.txt"},
    "tb_moon": {"ref": "03directTideMoon_icrf.txt"},
    "tb_sun": {"ref": "03directTideSun_icrf.txt"},
    "pole_tide": {"ref": "05poleTide_icrf.txt"},
    "ocean_pole_tide": {"ref": "06oceanPoleTide_icrf.txt"},
    "relativity": {"ref": "07relativistic_icrf.txt"},
    "dealiasing": {"ref": "08aod1b_RL06_icrf.txt"},
    "atmospheric_tide": {"ref": "09aod1b_atmosphericTides_icrf.txt"},
    "ocean_tide": {"ref": "11oceanTide_fes2014b_34major_icrf.txt"},
    "solid_earth_tide": {"ref": "04solidEarthTide_icrf.txt"},
}

if __name__ == "__main__":
    args = parser.parse_args()
    # eg:  -0.008047013077828 -0.008631299765786 0.018825345799512
    # tbs: 0.000000008974174 0.000000216144916 -0.000000176133478
    # tbm: -0.000000125185422 0.000000641617862 -0.000000440373072
    # rlt: -0.000000002813665 -0.000000003560216 0.000000015994246
    # sdt: -0.000000094344043 0.000000094020145 0.000000229112931
    # pt: 0.000000002147256 0.000000004491927 0.000000007623901
    # opt: 0.000000000382439 0.000000000726560 0.000000000740199
    # dls: 0.000000020517587 0.000000001453126 -0.000000025637531
    t, eg, tbs, tbm, rlt, sdt, pt, opt, dls, atm, ot = parse_acc()

    with PdfPages("foo.pdf") as pdf:
        fig1 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["earth_gravity"]["ref"]),
            t,
            eg,
            "m/sec^2",
            "Earth Gravity n,m>1",
        )
        pdf.savefig(fig1)
        plt.close(fig1)

        fig2 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["tb_moon"]["ref"]),
            t,
            tbm,
            "m/sec^2",
            "Third Body: Moon",
        )
        pdf.savefig(fig2)
        plt.close(fig2)

        fig3 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["tb_sun"]["ref"]),
            t,
            tbs,
            "m/sec^2",
            "Third Body: Sun",
        )
        pdf.savefig(fig3)
        plt.close(fig3)

        fig4 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["pole_tide"]["ref"]),
            t,
            pt,
            "m/sec^2",
            "(Solid Earth) Pole Tide",
        )
        pdf.savefig(fig4)
        plt.close(fig4)

        fig5 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["ocean_pole_tide"]["ref"]),
            t,
            opt,
            "m/sec^2",
            "Ocean Pole Tide, n,m>1",
        )
        pdf.savefig(fig5)
        plt.close(fig5)

        fig6 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["relativity"]["ref"]),
            t,
            rlt,
            "m/sec^2",
            "Relativistic Correction",
        )
        pdf.savefig(fig6)
        plt.close(fig6)

        fig7 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["dealiasing"]["ref"]),
            t,
            dls,
            "m/sec^2",
            "Dealiasing n,m>1",
        )
        pdf.savefig(fig7)
        plt.close(fig7)

        fig8 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["atmospheric_tide"]["ref"]),
            t,
            atm,
            "m/sec^2",
            "Atmospheric Tide",
        )
        pdf.savefig(fig8)
        plt.close(fig8)

        fig9 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["ocean_tide"]["ref"]),
            t,
            ot,
            "m/sec^2",
            "Ocean Tide",
        )
        pdf.savefig(fig9)
        plt.close(fig9)

        fig10 = genericPlotter(
            os.path.join(args.costg_dir, costg_data_files["solid_earth_tide"]["ref"]),
            t,
            sdt,
            "m/sec^2",
            "Solid Earth Tide",
        )
        pdf.savefig(fig10)
        plt.close(fig10)
