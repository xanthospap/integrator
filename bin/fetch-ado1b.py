#! /usr/bin/python

from ftplib import FTP

import datetime
import sys, os
import gzip
import shutil
import argparse

host = 'isdcftp.gfz.de'
gfz1Baod_dir = 'grace/Level-1B/GFZ/AOD'
valid_rls = {'RL05', 'RL06', 'RL07'}

def remoteFilename(rl, t, uncompressed=False):
    if rl.upper() not in valid_rls:
        raise RuntimeError(f'ERROR. Ivalid product number (RL) given: {rl}.')
    irl = int(rl[-1])
    ext = '.gz' if not uncompressed else ''
    return f"AOD1B_{t.strftime('%Y-%m-%d')}_X_{irl:02d}.asc{ext}"

def relDir(rl, t):
    if rl.upper() not in valid_rls:
        raise RuntimeError(f'ERROR. Ivalid product number (RL) given: {rl}.')
    return f"{rl.upper()}/{t.strftime('%Y')}"

def fetchFiles(ts, rl, local_dir, ftp_pass):
    fetched_files = []
    with FTP(host) as ftp:
        ftp.login(user='anonymous', passwd=ftp_pass)
    
        for t in ts:
            remote_fn = remoteFilename(rl, t)
            remote_dir = relDir(rl, t)
            target = os.path.join(os.path.join(gfz1Baod_dir,remote_dir), remote_fn)
            local_file = os.path.join(local_dir, remote_fn)

            with open(local_file, 'wb') as f:
                ftp.retrbinary(f'RETR {target}', f.write)

            if not os.path.isfile(local_file):
                raise RuntimeError(f"ERROR. Failed to download file {target}")
            
            fetched_files += [local_file]
    return fetched_files

def ugz(file):
    with gzip.open(file, 'rb') as f_in:
        with open(file[0:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def parseInDate(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d").date()
    except ValueError:
        raise argparse.ArgumentTypeError(f"Not a valid date: '{s}'. Expected format: YYYY-MM-DD")

class myFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description="",
    epilog=(
        """
    May, 2025"""
    ),
)

parser.add_argument(
    "-t",
    "--t-start",
    metavar="START_DATE",
    dest="t0",
    required=True,
    type=parseInDate,
    help="Start date of range of interest or date of interest if \'-e\' not specified.",
)

parser.add_argument(
    "-e",
    "--t-end",
    metavar="END_DATE",
    dest="t1",
    required=False,
    type=parseInDate,
    default=datetime.date.max,
    help="End date of interest (inclusive). If specified, we will download data for the range [START_DATE, END_DATE] (inclusive).",
)

parser.add_argument(
    "-r",
    "--product-type",
    metavar="RL",
    dest="rl",
    required=False,
    default="rl06",
    choices=valid_rls,
    help="The target product type, i.e. RL.",
)

parser.add_argument(
    "-p",
    "--ftp-pass",
    metavar="FTP_PASS",
    dest="ftp_pass",
    required=True,
    help="The ftp password. Note that we are using anonymous ftp so your email adress will do.",
)

parser.add_argument(
    "-d",
    "--dir",
    metavar="DOWNLOAD_DIR",
    dest="local_dir",
    required=False,
    default=os.getcwd(),
    help="The (local) dir where data shall be saved at.",
)

if __name__ == "__main__":
    args = parser.parse_args()

    delta = args.t1 - args.t0
    dates = [args.t0 + datetime.timedelta(days=i) for i in range(delta.days + 1)]

    skip_download_for = []
    for date in dates:
        if os.path.isfile(os.path.join(args.local_dir, remoteFilename(args.rl, date))) or os.path.isfile(os.path.join(args.local_dir, remoteFilename(args.rl, date, True))):
            skip_download_for += [date]
            print(f'Skipping download for file {remoteFilename(args.rl, date)}; file already present at {args.local_dir}')
            if not os.path.isfile(os.path.join(args.local_dir, remoteFilename(args.rl, date, True))):
                ugz(os.path.join(args.local_dir, remoteFilename(args.rl, date)))

    new_files = fetchFiles([t for t in dates if t not in skip_download_for], args.rl, args.local_dir, args.ftp_pass)

    for fl in new_files: ugz(fl)
