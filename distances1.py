"""
Script to flag and filter any duplicate breakpoints in same sample ID, or any adjacent breakpoints in same
sample ID.

Ruth McCole
November 7th 2016
"""


import argparse
import pandas as pd
import pybedtools as pbt



def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument('-f', "--file", type=str,
                        help='Filename for breakpoints file, with header, first four columns must be chr start stop ID')
    parser.add_argument('-d', '--distance', type=int, default=0,
                        help='Distance within which breakpoints are considered clustered, default is 0, so only '
                             'overlapping regions are clustered')

    return parser.parse_args()


def getDataIntoPandas(strFilename):
    #NB: this assumes your files DO have headers
    print 'Getting data from {0}'.format(strFilename)
    pdData = pd.read_csv(strFilename, sep='\t')
    return pdData


def pandaToBedtool(panda):
    arArFeatures = panda.values.tolist()
    btFeatures = getFeatures(arArFeatures)
    return btFeatures


def bedtoolToPanda(btobject):
    pdObject = pd.read_table(btobject.fn, header=None)
    return pdObject


def getFeatures(strFileName):
    btFeatures = pbt.BedTool(strFileName)
    return btFeatures


def sortBedLikePanda(pdData, arStrColumnNames):
    pdSorted = pdData.sort_values(by=arStrColumnNames)
    pdSorted.reset_index(inplace=True, drop=True)
    return pdSorted


def groupByCluster(pdClustered):
    groupCluster = pdClustered.groupby('cluster')

    arPdBreakpointsAloneInCluster = []

    arPdBreakpointsToKeepFromCluster = []

    arPdBreakpointsClusteredAndFilteredOut = []

    for tupName, pdGroup in groupCluster:
        if pdGroup.shape[0] == 1:
            #The group has only one breakpoint, keep it
            arPdBreakpointsAloneInCluster.append(pdGroup)

        if pdGroup.shape[0] > 1:
            #Check how many sampleIDs in group
            groupSample = pdGroup.groupby('sampleID')

            for tupName, pdSampleGroup in groupSample:
                if pdSampleGroup.shape[0] == 1:
                    #Only one row for this sample ID, keep it
                    arPdBreakpointsAloneInCluster.append(pdSampleGroup)
                elif pdSampleGroup.shape[0] > 1:
                    #There is more than one breakpoint from the same sample ID in this group. Keep the first one
                    arPdBreakpointsToKeepFromCluster.append(pdSampleGroup.iloc[[0], :])
                    #Also save the ones you are getting rid of
                    arPdBreakpointsClusteredAndFilteredOut.append(pdSampleGroup.iloc[1:])

    pdConcatBreakpointsAloneInCluster= pd.concat(arPdBreakpointsAloneInCluster, ignore_index=True)

    pdConcatBreakpointsToKeepFromCluster = pd.concat(arPdBreakpointsToKeepFromCluster, ignore_index=True)

    pdConcatBreakpointsClusteredAndFilteredOut = pd.concat(arPdBreakpointsClusteredAndFilteredOut, ignore_index=True)

    #Combine the two pandas that make up the ones you are keeping, and re-sort them.
    pdRowsToKeep = pd.concat([pdConcatBreakpointsAloneInCluster, pdConcatBreakpointsToKeepFromCluster], ignore_index=True)

    #Sort
    pdRowsToKeepSorted = sortBedLikePanda(pdRowsToKeep, ['chr', 'start'])

    #Count
    intTotalRows = pdClustered.shape[0]
    intRowsKept = pdRowsToKeepSorted.shape[0]
    intRowsFiltered = pdConcatBreakpointsClusteredAndFilteredOut.shape[0]

    if intTotalRows == intRowsKept + intRowsFiltered:
        print 'From {0} total rows, {1} rows were kept and {2} rows were filtered out, ' \
              'this adds up correctly'.format(intTotalRows, intRowsKept, intRowsFiltered)
    else:
        print "There is a PROBLEM. The rows don't add up:"
        print 'From {0} total rows, {1} rows were kept and {2} rows were filtered out, ' \
              'this adds up correctly'.format(intTotalRows, intRowsKept, intRowsFiltered)

    return pdRowsToKeepSorted, pdConcatBreakpointsClusteredAndFilteredOut


def savePanda(pdData, strFilename):
    pdData.to_csv(strFilename, sep='\t', index=False)


def main():
    args = get_args()

    intDistance = args.distance

    pdBreakpoints = getDataIntoPandas(args.file)

    btBreakpoints = pandaToBedtool(pdBreakpoints)

    #Sort
    btSorted = btBreakpoints.sort()

    #Produce clusters with bedtools cluster
    print 'Clustering breakpoints that are within {0} bp'.format(intDistance)
    btClustered = btSorted.cluster(d=intDistance)

    pdClustered = bedtoolToPanda(btClustered)

    pdClustered.columns = ['chr', 'start', 'end', 'sampleID', 'source', 'cluster']

    savePanda(pdClustered, 'Breakpoints_clustered_if_within_{0}bp.txt'.format(str(intDistance)))

    #For each cluster, split into clusters based on sample IDs, for each set of sample IDs with more than one member, keep just one
    print 'Filtering breakpoints'
    pdRowsToKeepSorted, pdConcatBreakpointsClusteredAndFilteredOut = groupByCluster(pdClustered)

    strClusteredFilename = '{0}_filtered_for_{1}_bp_clusters.txt'.format(args.file, intDistance)
    savePanda(pdRowsToKeepSorted, strClusteredFilename)

    strFilteredFilename = 'Rows_filtered_from_{0}_because_clustered_{1}_bp_away.txt'.format(args.file, intDistance)
    savePanda(pdConcatBreakpointsClusteredAndFilteredOut, strFilteredFilename)


if __name__ == "__main__":
    main()
