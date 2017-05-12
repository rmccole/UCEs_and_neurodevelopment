#!/usr/bin/env python

"""
Given two 3-column interval files, reports all bases overlapping between sets

Files should be formmatted as: chromosome/contig, start, stop seperated by tabs. Order does not matter, as only
reciprocal overlaps are reported. Make sure files are both 1-based starts and are collapsed. Overlaps may be missed
if intervals are not collapsed.

"""
import argparse


def get_args(strInput=None):
    parser = argparse.ArgumentParser("Reports all reciprocally overlapping bases between 2 given interval files. Files "
                                     "must be collapsed before use")
    parser.add_argument('A', type=argparse.FileType('rU'),
                        help="A 3-column interval file")
    parser.add_argument('B', type=argparse.FileType('rU'),
                        help="A 3-column interval file")
    if strInput:
        print "Given debug argument string: {0}".format(strInput)
        return parser.parse_args(strInput.split())
    return parser.parse_args()


def pointCheck(point, interval):
    """ Returns True if point is within interval

    Arguments:
    point = integer
    interval = chromosome/assemby, start, stop with start and stop being integers
    """
    assert isinstance(point, int)
    assert isinstance(interval[1], int)
    assert isinstance(interval[2], int)
    if interval[1] <= point <= interval[2]:
        return True
    return False


def formatInt(aInterval):
    """ Format an 3-column interval correctly """
    return [aInterval[0], int(aInterval[1]), int(aInterval[2])]


def intervalCheck(intervalA, intervalB):
    """Return an interval containing all bases in common between intervalA and intervalB

    Takes an interval in 3-column format: chr, start, stop
    """
    chrA, startA, stopA = intervalA
    chrB, startB, stopB = intervalB
    if chrA == chrB:
        overlapChr = chrA
        # Check if start coordinate of interval A lies within interval B
        if pointCheck(startA, intervalB):
            overlapStart = startA
            if stopA <= stopB:
                overlapStop = stopA
            else:
                overlapStop = stopB
            return [overlapChr, overlapStart, overlapStop]
            # If not, check if end coordinate of interval A lies within interval B
        elif pointCheck(stopA, intervalB):
            overlapStop = stopA
            overlapStart = startB
            return [overlapChr, overlapStart, overlapStop]
        # If not, check if interval A surrounds interval B
        elif startA < startB and stopA > stopB:
            overlapStart, overlapStop = startB, stopB  # Report smaller of the two
            return [overlapChr, overlapStart, overlapStop]
        else:
            return False
    return False


def main(args):
    fileA = [formatInt(line.strip().split('\t')) for line in args.A]
    fileA.sort(key=lambda x: (x[0], x[1], x[2]))
    fileB = [formatInt(line.strip().split('\t')) for line in args.B]
    fileB.sort(key=lambda x: (x[0], x[1], x[2]))
    overlaps = []
    for intervalB in fileB:
        for intervalA in fileA:
            if intervalB[0] != intervalA[0]:  # Check chromosome is the same for both intervals
                continue
            elif intervalB[1] > intervalA[2]:  # Don't check against intervalA if stop is less than intervalB start
                continue
            overlap = intervalCheck(intervalB, intervalA)
            if overlap:  # If overlap found, append and check next intervalB
                overlaps.append(overlap)
                continue
            elif intervalB[2] < intervalA[1]:  # Don't check against intervalA if intervalB ends before intervalA starts
                break
            else:
                break
    return overlaps


if __name__ == '__main__':
    args = get_args()
    overlaps = main(args)
    for line in overlaps:
        print "\t".join(map(str, line))