#!/usr/bin/env python

import pysam
import argparse
import sys
import time
import collections as c
import math
import socket
import re
from multiprocessing import Pool
import resource

if socket.gethostname() == 'SEQ-SRV-03':
    print("Working on SEQ-SRV-03 - fixing path...")
    for n, i in enumerate(sys.path):
        if re.match(".*pandas.*", i) is not None:
            del(sys.path[n])

import pandas as pd
import numpy as np


__author__ = 'kdavie'

""" Takes a sam or bam file as an input and outputs a matrix of either
read depth or fragment start sites based off of a bed input. An aggregation
plot and a joined aggregation plot/heatmap will be produced in two separate
files alongside the matrix. Multithreading is enabled drastically reducing
the time required to generate the matrix. Plot generation takes a small
amount of time.
"""

# Help and arguments parser
parser = argparse.ArgumentParser(description='Takes a sam or bam file as an input and outputs a 6 column bed file '
                                             'with the locations of each read start point and its strand, col 5 = . '
                                             'Also able to create a matrix of either read depth or fragment start sites'
                                             ' based off of a bed input. When the matrix option is used, a matrix will '
                                             'be created based on the defined parameters, this can easily be used to '
                                             'create heatmaps and aggregation plots.')
parser.add_argument('inputFile', metavar='input', type=str, help='Input data file - Either SAM, BAM or Tabix formatted '
                                                                 'bed file (autodetected)', default=sys.stdin)
parser.add_argument('bed', metavar='bed', type=str, help='Bed file for creating a matrix')
parser.add_argument('size', metavar='size', type=int, help='How wide should we make the matrix from each peak?',
                    default=0)
parser.add_argument('outputFile', metavar='output', type=str, help='Output file - either bed file or matrix',
                    default=sys.stdout)
parser.add_argument('--quiet', action='store_true', help='Outputs no information || false', default=False)
parser.add_argument('--atac', action='store_true', help='Processes sites with a 9bp difference to account for binding'
                                                        ' vs cutting in ATAC-seq || false', default=False)
parser.add_argument('-pos', metavar='bp', type=int, help='Override default ATAC-seq positive shift || 4',
                    default=0)
parser.add_argument('-neg', metavar='bp', type=int, help='Override default ATAC-seq Negative shift || -5',
                    default=0)
parser.add_argument('--debug', action='store_true', help='Print debug info, confirm every read', default=False)
parser.add_argument('-smooth', metavar='bp', type=int, help='Smoothening to apply in bp, only for cutsites', default=1)
parser.add_argument('-binbed', action='store_true', help='Process bed as binary?')
parser.add_argument('--sumsOnly', action='store_true', help='Only provide the sums of the matrix')
parser.add_argument('--rpm', action='store_true', help='Normalize the matrix to reads per million')
parser.add_argument('--regNorm', action='store_true', help='Normalize the matrix to number of regions')
parser.add_argument('-ylim', type=float, help='Set the upper y limit (lower will be 0)')
parser.add_argument('-proc', type=int, help='Set the number of processes to use || 1', default=1)
parser.add_argument('-quantile', type=float, help='Set the quantile to normalise by || 0.95', default=0.95)
parser.add_argument('--singleCell', action='store_true', help='Not implemented', default=False)

args = parser.parse_args()

if args.debug:
    print("\n\nDebug info:\n")
    print(" Input file:             " + args.inputFile)
    print(" Output file:            " + args.outputFile)
    print(" Bed file used           " + str(args.bed))
    print(" Quiet?:                 " + str(args.quiet))
    print(" ATAC-Seq?               " + str(args.atac))
    print(" Manual Positive shift   " + str(args.pos))
    print(" Manual Negative shift   " + str(args.neg))
    print(" Smoothing window        " + str(args.smooth))
    print(" Bed vs Bed as binary?   " + str(args.binbed) + '\n\n')

# Try and load the data file - Can work with both B/SAM and Tabix format
# TODO: Add support for bigWig files - Any other we might need?

try:
    datafile = pysam.Samfile(args.inputFile, "rb")
    isSam = True
    isTabix = False

# TODO: Emtpy sam files appear as tabix files
except ValueError:
    try:
        datafile = pysam.Tabixfile(args.inputFile, "r")
        isTabix = True
        isSam = False
        if not args.quiet:
            print("Input is bed format")
    except ValueError:
        print("Input file is of an unsupported type")
        exit(1)
except OSError:
    print("ERROR: Input file not found! Exiting.")
    exit(1)

# Try and load the bed file

try:
    b = open(args.bed, 'r')
except TypeError:
    print("Cannot open bed file")
    exit(1)
except OSError:
    print("ERROR: Bed file not found! Exiting")
    exit(1)

# Assign some initial variables
outBed = []
reverse = 0
forward = 0
percent = 0
bufferPos = 0
bufferNeg = 0
bedLen = 0
zerr = 0
mSums = ["Sums"]

if args.rpm and not isSam:
    print("Cannot RPM normalize a bed file - will continue without normalisation")
    args.rpm = False
if args.rpm and args.atac:
    print("Cannot RPM normalize a cutsites - will continue without normalisation")
    args.rpm = False
if args.rpm:
    libSize = datafile.mapped

# Deal with ATAC overrides
if args.atac:
    if args.pos > 0:
        bufferPos = args.pos
    else:
        bufferPos = 4
    if args.neg > 0:
        bufferNeg = args.neg
    else:
        bufferNeg = -5


# Function for making a matrix from read depth
def makematrix(region):

    try:
        datafile = pysam.Samfile(args.inputFile, "rb")
        isTabix = False
        isSam = True
    except ValueError:
        datafile = pysam.Tabixfile(args.inputFile, "r")
        isTabix = True
        isSam = False

    count = 0

    # Process the region name and set some initial variables

    if args.debug:
        print("Working on - " + region)

    chrom = region.split('\t')[0]
    start = int(region.split('\t')[1])
    end = int(region.split('\t')[2])
    count += 1
    cnt = c.defaultdict(lambda: 0)
    totalreads = 0
    lineName = chrom + ":" + str(start) + "-" + str(end)

    # Change the start and end sites to be Half way between +/- size
    oldstart = start
    oldend = end

    start = int((oldstart + (math.floor((oldend - oldstart) / 2))) - args.size)
    end = int((oldstart + (math.floor((oldend - oldstart) / 2))) + args.size)
    matrixline = []

    # Try and extract strand information, if non existent or not '+'/'-', assume positive strand
    try:
        strand = region.split('\t')[5].rstrip('\n')
        if strand != '+' or strand != '-':
            strand = '+'
    except IndexError:
        strand = '+'

    startpileup = time.time()

    # Support for Tabix indexed bed files
    if isTabix:
        try:
            for reg in datafile.fetch(chrom, start, end):
                sreg = reg.split('\t')
                if sreg[0] == chrom and (start <= int(sreg[1]) <= end or
                                         start <= int(sreg[2]) <= end):
                    for i in range(int(sreg[1]), int(sreg[2]), 1):
                        if args.binbed:
                            cnt[i] += 1
                        try:
                            cnt[i] += int(float(sreg[4]))
                        except (IndexError, ValueError):
                            cnt[i] += 1
        except ValueError:
            if args.debug:
                print("Region " + chrom + " " + str(start) + " " + str(end) + " does not exist in this file")
            datafile.close()
            matrixline = pd.Series(np.zeros((args.size * 2) + 1), name=lineName)
            return (matrixline)

        if args.debug:
            print("Finished checking bed in " + str(time.time() - startpileup) + " seconds")

        if strand == '+':
            for base in range(start, end + 1, 1):
                if base < start or base > end:
                    continue
                totalreads += int(cnt[base])
                matrixline.append(float(cnt[base]).rstrip('\n'))

        if strand == '-':
            for base in range(end + 1, start, 1):
                if base < start or base > end:
                    continue
                totalreads += int(cnt[base])
                matrixline.append(float(cnt[base]).rstrip('\n'))
            matrixline = matrixline[::-1]

    # Support for SAM and BAM files
    if isSam:
        cols = 0
        nextcol = start

        if start < 0:
            start = 0

        for pileupColumn in datafile.pileup(chrom, start, end):
            if int(start) > int(pileupColumn.pos) or int(end) < int(pileupColumn.pos):
                continue
            while nextcol < pileupColumn.pos:
                matrixline.append(0)
                nextcol += 1
                cols += 1
            matrixline.append(float(pileupColumn.n))
            if args.debug:
                print("Appending " + str(float(pileupColumn.n)))
            totalreads += pileupColumn.n
            nextcol += 1
            cols += 1

        if args.debug:
            print("Total number of pileup columns: " + str(cols))
            print("Command sent: reads.pileup(" + str(chrom) + ", " + str(start) + ", " + str(end) + ")")
            print("cols - (end - start) = " + str(cols - (end - start)))
            print(str(nextcol) + ", " + str(end) + ", " + str(start) + ", " + str(cols))

        if nextcol <= end:
            for i in range(nextcol, end + 1, 1):
                matrixline.append(0)

        if strand == '-':
            matrixline = matrixline[::-1]

        datafile.close()

        if len(matrixline) == ((args.size * 2) + 1):
            matrixline = pd.Series(matrixline, name=lineName)
            return(matrixline)
        else:
            matrixline = pd.Series(np.zeros((args.size * 2) + 1), name=lineName)
            return(matrixline)


def cutsiteMatrix(region):

    try:
        datafile = pysam.Samfile(args.inputFile, "rb")
        isTabix = False
        isSam = True
    except ValueError:
        datafile = pysam.Tabixfile(args.inputFile, "r")
        isTabix = True
        isSam = False

    count = 0
    # Split the input bed
    if args.debug:
        print("Working on - " + region)
    chrom = region.split('\t')[0]
    start = int(region.split('\t')[1])
    end = int(region.split('\t')[2])
    count += 1
    cnt = c.defaultdict(lambda: 0)
    lineName = chrom + ":" + str(start) + "-" + str(end)

    # If a size is specified change the start and end sites to be Half way between +/- size, if not
    oldstart = start
    oldend = end

    start = (int((oldstart + (math.floor((oldend - oldstart) / 2))) - args.size) - args.smooth)
    end = (int((oldstart + (math.floor((oldend - oldstart) / 2))) + args.size) + args.smooth)
    matrixline = []

    # Try and extract strand information, if non existent or not '+'/'-', assume positive strand
    try:
        strand = region.split('\t')[5].rstrip('\n')
        if strand != '+' or strand != '-':
            strand = '+'
    except IndexError:
        strand = '+'

    startpileup = time.time()
    cellsSeen = set()
    # Support for Tabix indexed bed files
    if isTabix:
        print('Tabix cutsites not yet supported')
        exit(1)

        #TODO: Reimplement this
        try:
            for reg in datafile.fetch(chrom, start, end):
                sreg = reg.split('\t')
                if sreg[0] == chrom and (start <= int(sreg[1]) <= end or
                                         start <= int(sreg[2]) <= end):
                    cnt[sreg[1]] += 1
        except ValueError:
            if args.debug:
                print("Region " + chrom + " " + str(start) + " " + str(end) + " does not exist in this file")
            datafile.close()
            matrixline = pd.Series(np.zeros((args.size * 2) + 1), name=lineName)
            return(matrixline)

        if args.debug:
            print("Finished checking bed in " + str(time.time() - startpileup) + " seconds")

        if strand == '+':
            for base in range(start, end + 1, 1):
                if base < start or base > end:
                    continue
                totalreads += int(cnt[base])
                matrixline.append(float(cnt[base]).rstrip('\n'))

        if strand == '-':
            for base in range(end + 1, start, 1):
                if base < start or base > end:
                    continue
                totalreads += int(cnt[base])
                matrixline.append(float(cnt[base]).rstrip('\n'))
            matrixline = matrixline[::-1]

    # Support for SAM and BAM files
    if isSam:
        if start < 0:
            start = 0

        for read in datafile.fetch(chrom, start, end):
            if args.singleCell:
                cellName = read.query_name.split(':')[0]
                if cellName in cellsSeen:
                    break
                cellsSeen.add(cellName)
            if not read.is_reverse:
                if start <= int(read.pos) <= end:
                    cnt[read.pos] += 1
            if read.is_reverse:
                if start <= int(read.pos + read.inferred_length) <= end:
                    cnt[read.pos + read.inferred_length] += 1
        for bp in range(start + 1, end):
            matrixline.append(cnt[bp])

        datafile.close()
        if len(matrixline) == ((args.size * 2) + 1):
            matrixline = pd.Series(matrixline, name=lineName)
            return(matrixline)
        else:
            matrixline = pd.Series(np.zeros((args.size * 2) + 1), name=lineName)
            return(matrixline)

if not args.quiet:
    print("\nBeginning to process file\n")

index = []
bedLines = []

# Check to see if any of the chromosomes start with chr. If they do, we assume
# they all do, we will fix the bed file to be the same.

chrRef = False

for ref in datafile.references:
    if ref.startswith('chr'):
        chrRef = True

# Check and process the regions in the bed file

for line in b:
    line = line.rstrip('\n').split('\t')
    [str(x) for x in line]

    # Fix the bed regions to be the same reference names as the datafile
    if line[0].startswith('chr') and not chrRef:
        line[0] = line[3:]
    elif not line[0].startswith('chr') and chrRef:
        line[0] = 'chr' + line[0]

    lineName = line[0] + ':' + line[1] + '-' + line[2]

    if len(line) >= 3:
        if lineName not in index:
            bedLen += 1
            bedLines.append("\t".join(line))
            index.append(lineName)
        else:
            if not args.quiet:
                print("Duplicate entry found - " + "\t".join(line))
    elif len(line) >= 1:
            if not args.quiet:
                print("Malformed line - " + "\t".join(line))

cols = [x for x in range(-args.size, args.size + 1)]

mlines = []
start_time = time.time()

# Start processing all regions across multiple cores if needed
if args.debug:
    print('4. Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


p = Pool(args.proc)

if args.atac:
    dicts = p.map(cutsiteMatrix, bedLines)
else:
    dicts = p.map(makematrix, bedLines)

p.close()
p.join()

if not args.quiet:
    print("Generating matrix...")
# matrix = pd.DataFrame.from_dict(allData, orient='index')
matrix = pd.concat(dicts, axis=1)
del(dicts)
matrix = matrix.fillna(0)
matrix = matrix.astype(float)
matrix = matrix.T

# del(allData)
matrix.columns = cols
matrix = matrix.reindex(index)
del(index)

if args.debug:
    print('5. Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


if args.rpm:
    if not args.quiet:
        print("Normalising matrix...")
    matrix = (matrix.fillna(0) / libSize) * 1000000

if args.regNorm:
    matrix = matrix / int(bedLen)

if args.smooth > 1:
    if not args.quiet:
        print("Smoothing matrix...")
    matrix = pd.rolling_mean(matrix.T, args.smooth).T
    if not args.quiet:
        print("Selecting columns")
    smoothMatrix = pd.DataFrame(index=index)
    for i in matrix.columns:
        if args.debug:
            print("On column" + str(i))
        if int(i) % args.smooth == 0 or int(i) == 0:
            if int(i) == args.size or int(i) == -(args.size):
                continue
            smoothMatrix[i] = matrix[i]
    matrix = smoothMatrix

if not args.quiet:
    print("Calculating sums...")

if args.debug:
    print('6. Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

matrixSums = pd.DataFrame(matrix.sum()).T

if not args.quiet:
    print("Generating plot...")


x_ticks = [-args.size, int(-(args.size/2)), 0, int((args.size/2)), args.size]
x = [i + args.size for i in x_ticks]
plt.figure(figsize=(5, 15))
plt.subplot2grid((10, 1), (0, 0))
plt.plot(matrix.sum(), c='darkorange')
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.xticks()
plt.grid(True)
plt.subplot2grid((10, 1), (0, 1), rowspan=9)
plt.grid(False)
plt.xticks(x, x_ticks)
if not args.atac:
    plt.imshow(matrix,
               interpolation='nearest',
               vmax=matrix.quantile(q=args.quantile, axis=1).mean(),
               cmap='Oranges',
               aspect='auto')
else:
    plt.imshow(matrix,
               interpolation='nearest',
               vmax=1,
               cmap='Oranges',
               aspect='auto')

plt.savefig('.'.join([args.outputFile, 'heatmap', 'pdf']), bbox_inches='tight')


plt.figure()
mpl.use('Agg')
import matplotlib.pyplot as plt
if args.smooth > 1:
    matrix.sum().plot(xlim=(-(int(args.size)) + 10, (int(args.size)) - 10), color="black")
else:
    matrix.sum().plot(xlim=(-(int(args.size)), (int(args.size))), color="black")
if args.rpm:
    plt.ylabel('Aggregate reads per (mapped) million')
else:
    plt.ylabel('Aggregate reads')

plt.xlabel('Position')

if args.ylim is not None:
    plt.ylim(0, args.ylim)
plt.axvline(x=0, ls='--', lw=1)
plt.savefig('.'.join([args.outputFile, 'aggPlot', 'pdf']), bbox_inches='tight')

if args.rpm:
    matrixSums.index = ['Aggregate Reads per (mapped) million']
else:
    matrixSums.index = ['Aggregate Reads']
if args.sumsOnly:
    matrixSums.to_csv(args.outputFile, sep='\t', na_rep='0')
else:
    if not args.quiet:
        print("Writing full matrix. This can take some time and a lot of space.")
    matrix.to_csv(args.outputFile, sep='\t', na_rep='0')
