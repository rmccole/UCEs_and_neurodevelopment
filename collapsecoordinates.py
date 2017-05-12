#!/usr/bin/env python
"""
Collapses overlapping coordinates of a given interval file. Sorts file first by
chr, start and stop. Then, runs through each interval, storing the previous one
and taking the next one in. It then compares the two for any overlap, adjusts
the interval and prints it.

Chamith Fonseka
7 February 2011
"""
import argparse
import sys


def get_args(strInput=False):
    parser = argparse.ArgumentParser(description="Collapses overlapping " +
                                                 "coordinates of a 1-based interval file.")
    parser.add_argument("file", type=argparse.FileType('rU'),
                        help="A 3-column interval file")
    if strInput:
        return parser.parse_args(strInput.split())
    else:
        return parser.parse_args()


def stdout_writer(aList):
    for line in aList:
        print '\t'.join(map(str, line))


def collapse(aIntervals):
    # Initialize variables
    strChr = iStart = iStop = 0
    aOut = []
    # Write status message
    sys.stderr.write("Collapsing coordinates...\n")
    for aInterval in aIntervals:
        # Test if an interval has been stored (always past first loop)
        if strChr:
            # Test if next interval is on a different chr OR if start of
            # next interval is larger than stop of previous interval
            if strChr != aInterval[0] or aInterval[1] > (iStop + 1):
                # Write interval and clear start/stop coordinates
                aOut.append([strChr, iStart, iStop])
                iStart = iStop = 0
        strChr = aInterval[0]
        # Advance to next interval if iStart is empty
        if not iStart:
            iStart = aInterval[1]
            # If next interval overlaps, adjust stop to larger coordinate
        if aInterval[2] > iStop:
            iStop = aInterval[2]
            # Write last line
    aOut.append([strChr, iStart, iStop])
    return aOut


if __name__ == "__main__":
    args = get_args()
    aIntervals = map(lambda x: [x[0], int(x[1]), int(x[2])],
                     [line.strip().split('\t') for line in args.file])
    aIntervals.sort(key=lambda x: (x[0], x[1], x[2]))
    aOut = collapse(aIntervals)
    stdout_writer(aOut)
