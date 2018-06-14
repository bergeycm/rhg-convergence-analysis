#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Fake IMPUTE parts, pretending they were created just before full IMPUTE results
# --- This allows us to skip their transfer with rsync and just focus on the full ones.
# ----------------------------------------------------------------------------------------

import re
import os

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

pattern = re.compile("^chr([0-9]+)$")

orig_times = []

# Find earliest create IMPUTE file
with open("data/hg19.genome") as gen:
    for line in gen:
        chr = line.split()[0]
        m = pattern.match(chr)
        if (m):
            chr_num = m.group(1)
            chr_len = int(line.split()[1])
            orig = "results/impute/AGRHUM_EASTERN_100x267251.1M." + str(chr_num) + ".1000g.gen.impute2"
            orig_times.append(os.stat(orig).st_mtime)

earliest = min(orig_times)

with open("data/hg19.genome") as gen:
    for line in gen:
        chr = line.split()[0]
        m = pattern.match(chr)
        if (m):
            chr_num = m.group(1)
            chr_len = int(line.split()[1])
            orig = "results/impute/AGRHUM_EASTERN_100x267251.1M." + str(chr_num) + ".1000g.gen.impute2"
            prefix = orig + "_pt"
            for i in range(0, 1 + chr_len // 5000000):
                pt_name = prefix + str(i)
                touch(pt_name, (earliest - 60, earliest - 60))
