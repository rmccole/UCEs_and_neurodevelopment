"""
Script to analyze two files of genome features and get information about their relative distances.


Ruth McCole
September 7th 2016
"""

import argparse
import pandas as pd
import pybedtools as pbt
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument('-q', "--queryFile", type=str,
                        help='Filename for first feature file, with header, should be ZERO based, these are to be matched'
                             'with the random files, if using')
    parser.add_argument('-ql', '--queryLabel', type=str,
                        help='Label for query')
    parser.add_argument('-r', '--refFile', type=str,
                        help='Filename for second feature file, with header, should be ZERO based')
    parser.add_argument('-rl', '--refLabel', type=str,
                        help='Label for reference')
    parser.add_argument('-f', '--filenames', type=str,
                        help='File containing filenames of random files, these are in place of query file, and are ONE BASED')
    parser.add_argument('-b', '--bins', type=int, default=50,
                        help='Number of bins to aggregate into')
    parser.add_argument('-x', '--xLims', type=int, nargs=2, default=(0, 20000),
                        help='Upper and lower limits of bins for data aggrgation, supply two integers, separated by'
                             'a space')
    parser.add_argument('-k', '--keepTies', action='store_true', dest='boolK',
                        help='pass -k if you would like to keep all instances where a single query region is '
                             'equidistant from two reference items. This can happen for example when a UCE'
                             'is in a region that is annotated as 2 different genes.')

    return parser.parse_args()


def getFeatures(strFileName):
    btFeatures = pbt.BedTool(strFileName)
    return btFeatures


def getDataIntoPandas(strFilename):
    #NB: this assumes your files DO have headers
    print 'Getting data from {0}'.format(strFilename)
    pdData = pd.read_csv(strFilename, sep='\t')
    return pdData


def getDataIntoPandasNoHeader(strFilename, arHeaders):
    #This is for files without headers
    #print 'Getting data from {0}'.format(strFilename)
    pdData = pd.read_csv(strFilename, sep='\t', header=None, names=arHeaders)
    return pdData


def get_Filenames_file(strFilenames):
    """
    Obtain a list of the filenames of the files to process
    """
    with open(strFilenames, 'r') as Filenames_filehandle:
        aFiles = [line.strip() for line in Filenames_filehandle]
    return aFiles


def getPdRandom(strRandomFilename):
    pdRandom = getDataIntoPandasNoHeader(strRandomFilename, ['chr', 'start', 'end'])
    return pdRandom


def makeRandom0based(pdRandom):
    #Make 0 based
    pdRandom['0based_start'] = pdRandom.apply(lambda row: (row['start'] - 1), axis=1)
    pdRandom0based = pdRandom[['chr', '0based_start', 'end']]

    pdRandom0based.columns = ['chr', 'start', 'end']

    return pdRandom0based


def get_randoms_make_0based(arFilenames):
    arPdRandoms = []
    for strFilename in arFilenames:
        pdRandom = getPdRandom(strFilename)
        pdRandom0based = makeRandom0based(pdRandom)
        arPdRandoms.append(pdRandom0based)

    return arPdRandoms


def pandaToBedtool(panda):
    arArFeatures = panda.values.tolist()
    btFeatures = getFeatures(arArFeatures)
    return btFeatures


def bedtoolToPanda(btobject):
    pdObject = pd.read_table(btobject.fn, header=None)
    return pdObject


def savePanda(pdData, strFilename):
    pdData.to_csv(strFilename, sep='\t', index=False)


def saveBedTool(btObject, strFilename):
    btObject.saveas(strFilename)


def getChrSetFromPanda(pdData):
    #Assumes the chr column is the first column
    arChrs = pdData.iloc[:, 0].values.tolist()
    setChrs = set(arChrs)
    return setChrs


def renameFirstColumnToChr(pdData):
    #List current column headers
    arPdColumns = list(pdData)
    #Find how many there are
    intLenPdColumns = len(arPdColumns)
    #Create a new list of headers, replacing the first with 'chr', and using the rest as is
    arNewColumns = ['chr'] + arPdColumns[1:intLenPdColumns]
    #Assign these to the dataframe
    pdData.columns = arNewColumns
    return pdData


def filterPdUnmatchedChr(pdBreakpoints, pdUCEs):
    #Want to filter and save any breakpoints or UCEs that are on a chromosome without a
    #corresponding feature in the other file.

    #Make set of chromosomes
    setBreakpointChrs = getChrSetFromPanda(pdBreakpoints)
    setUCEChrs = getChrSetFromPanda(pdUCEs)

    #Find breakpoint chrs that are not in UCE chrs
    setOnlyBreakpointChrs = setBreakpointChrs.difference(setUCEChrs)
    setOnlyUCEChrs = setUCEChrs.difference(setBreakpointChrs)
    #print 'The following chromosomes are in dataset A but not B {0}'.format(str(setOnlyBreakpointChrs))
    #print 'The following chromosomes are in dataset B but not A {0}'.format(str(setOnlyUCEChrs))

    #Now find all entries in data that correspond to these chromosomes, save them, and filter them from the data
    #going forward
    pdBreakpointsColumnNamed = renameFirstColumnToChr(pdBreakpoints)

    #Get only the rows where the entry in the first column (called 'chr' matches the only breakpoint chromosomes
    pdBreakpointsUnmatchedRows = pdBreakpointsColumnNamed[pdBreakpointsColumnNamed.chr.isin(setOnlyBreakpointChrs)]
    #save
    savePanda(pdBreakpointsUnmatchedRows, 'Breakpoints_without_a_UCE_on_that_chr.txt')

    #Find UCE chrs that are not in breakpoint file
    pdUCEsColumnNamed = renameFirstColumnToChr(pdUCEs)

    pdUCEsUnmatchedChrs = pdUCEsColumnNamed[pdUCEsColumnNamed.chr.isin(setOnlyUCEChrs)]
    savePanda(pdUCEsUnmatchedChrs, 'UCEs_without_a_breakpoint_on_that_chr.txt')

    #Now filter the data to get rid of these rows
    pdBreakpointsWithUCEs = pdBreakpointsColumnNamed[~pdBreakpointsColumnNamed.chr.isin(setOnlyBreakpointChrs)]

    pdUCEsWithBreakpoints = pdUCEsColumnNamed[~pdUCEsColumnNamed.chr.isin(setOnlyUCEChrs)]

    return pdBreakpointsWithUCEs, pdUCEsWithBreakpoints


def getClosestAndSave(btFileredBreakpoints, btFilteredUCEs):
    #First sort
    btFilteredBreakpointsSorted = btFileredBreakpoints.sort()
    btFilteredUCEsSorted = btFilteredUCEs.sort()

    #Get closest UCEs to breakpoints
    btBreakpointsWithClosestUCE = btFilteredBreakpointsSorted.closest(btFilteredUCEsSorted, d=True)

    #Get closest breakpoints to UCEs
    btUCEsWithClosestBreakpoints = btFilteredUCEsSorted.closest(btFilteredBreakpointsSorted, d=True)

    return btBreakpointsWithClosestUCE, btUCEsWithClosestBreakpoints


def getSortedBedtool(strFilename):
    pdFeatures = getDataIntoPandas(strFilename)
    btFeatures = pandaToBedtool(pdFeatures)
    btSortedFeatures = btFeatures.sort()
    return btSortedFeatures


def AddHeadersSortClosestResults(btClosestResult, arHeaders):
    pdClosestBreakpoints = bedtoolToPanda(btClosestResult)
    pdClosestBreakpoints.columns = arHeaders

    pdClosestBreakpointsSortedDistances = pdClosestBreakpoints.sort(columns='Distance')
    savePanda(pdClosestBreakpointsSortedDistances, 'Breakpoints_with_closest_UCEs.bed')

    return pdClosestBreakpointsSortedDistances


def AddKbMbColumns(pdDistances):
    #Obtain the distances scaled to kb
    pdDistances['Distance_Kb'] = pdDistances.apply(lambda row: (float(row['Distance'])/1000), axis=1)
    #Obtain the distances in Mb
    pdDistances['Distance_Mb'] = pdDistances.apply(lambda row: (float(row['Distance'])/1000000), axis=1)

    return pdDistances


def prepareClosestResults(btClosestResult, arHeaders):
    #First, add headers and sort
    pdClosestSorted = AddHeadersSortClosestResults(btClosestResult, arHeaders)
    #Then, add columns with different distance units
    pdWithDistanceColumns = AddKbMbColumns(pdClosestSorted)
    #Then, remove any rows with Distance = -1
    #First create a mask
    pdDistanceNotMinusOne = pdWithDistanceColumns.Distance != -1
    #Apply the mask
    pdBreakpointsSortedDistancesNoMinusOnes = pdWithDistanceColumns[pdDistanceNotMinusOne]

    return pdBreakpointsSortedDistancesNoMinusOnes


def drawAndSaveHistogram(pdData, strColumn, strFilename, strXLabel):
    seriesDistancesMb = pdData[strColumn]
    plt.figure()
    sns.set_style('whitegrid')
    g = sns.distplot(seriesDistancesMb, kde=False, rug=False)
    sns.axlabel(strXLabel, 'Counts')
    plt.savefig(strFilename, bbox_inches='tight')
    plt.clf()


def GetFeaturesWithinCertainDistance(pdBreakpointsWithClosestUCESorted, intDistanceKb):
    #Int to float
    floatDistanceKb = float(intDistanceKb)
    #Get boolean mask
    boolDistanceLessThanIntDistance = pdBreakpointsWithClosestUCESorted.Distance_Kb < floatDistanceKb
    #Apply mask
    pdOnlyUnderSpecifiedDistance= pdBreakpointsWithClosestUCESorted[boolDistanceLessThanIntDistance]
    #Get number of elements in dataframe
    intNumberBreakpointsWithUCEsUnderDistance = len(pdOnlyUnderSpecifiedDistance)

    return intNumberBreakpointsWithUCEsUnderDistance, pdOnlyUnderSpecifiedDistance


def zoomedHistogram(pdData, strColumn, strFilename, strXLabel, intXLim, bins):
    seriesDistancesMb = pdData[strColumn]
    plt.figure()
    sns.set_style('whitegrid')
    g = sns.distplot(seriesDistancesMb, kde=False, rug=False, bins=bins)
    sns.plt.xlim(0, intXLim)
    sns.axlabel(strXLabel, 'Counts')
    plt.savefig(strFilename, bbox_inches='tight')
    plt.clf()


def convertToPandaAddHeaderSortOnDistance(btData, arHeaders):
    pdData = bedtoolToPanda(btData)
    pdData.columns = arHeaders

    #Sort by distance column
    pdDataSorted = pdData.sort_values(by='Distance')
    return pdDataSorted


def generateHeaders(pdData, strStem):
    intLenQuery = pdData.shape[1]
    arFirst3Headers = ['chr_' + strStem, 'start_' + strStem, 'end_' + strStem]
    arNextHeaders = []
    for integer in range(3, intLenQuery):
        nextHeader = strStem + str(integer)
        arNextHeaders.append(nextHeader)
    arAllHeaders = arFirst3Headers + arNextHeaders

    return arAllHeaders


def mainGetDistances(pdBreakpoints, pdUCEs, boolKeepQueryTies):
    #Save headers for later
    arFilteredBreakpointHeaders = list(pdBreakpoints)
    arFilteredUCEHeaders = list(pdUCEs)

    #Filter for unmatched chromosomes, and change first column to 'chr'
    pdFilteredBreakpoints, pdFilteredUCEs = filterPdUnmatchedChr(pdBreakpoints, pdUCEs)

    #Get into bedtools
    btFilteredBreakpoints = pandaToBedtool(pdFilteredBreakpoints)
    btFilteredUCEs = pandaToBedtool(pdFilteredUCEs)

    #Find and save closest elements
    btBreakpointsWithClosestUCE, btUCEsWithClosestBreakpoints = getClosestAndSave(btFilteredBreakpoints, btFilteredUCEs)

    #Prepare headers
    ###Version 4 - provide headers so that no two are the same when query and refs are combined
    arBreakpointsGeneratedHeaders = generateHeaders(pdFilteredBreakpoints, 'Ref')
    arQueryGeneratedHeaders = generateHeaders(pdFilteredUCEs, 'Query')

    arNewBreakpointHeaders = arBreakpointsGeneratedHeaders + arQueryGeneratedHeaders + ['Distance']
    arNewUCEHeaders = arQueryGeneratedHeaders + arBreakpointsGeneratedHeaders + ['Distance']

    #arNewBreakpointHeaders = arFilteredBreakpointHeaders + arFilteredUCEHeaders + ['Distance']

    #arNewUCEHeaders = arFilteredUCEHeaders + arFilteredBreakpointHeaders + ['Distance']

    #Convert to pandas with these headers, and sort by distance
    pdBreakpointsSortedOnDistance = convertToPandaAddHeaderSortOnDistance(btBreakpointsWithClosestUCE, arNewBreakpointHeaders)
    pdUCEsSortedOnDistance = convertToPandaAddHeaderSortOnDistance(btUCEsWithClosestBreakpoints, arNewUCEHeaders)

    #Add distance columns with other units
    pdBreakpointsWithDistanceUnits = AddKbMbColumns(pdBreakpointsSortedOnDistance)
    pdUCEsWithDistanceUnits = AddKbMbColumns(pdUCEsSortedOnDistance)

    ###Version 3, remove full duplicates, for query file only. NOT DONE FOR REFS.
    pdUCEsWithDistanceUnits.drop_duplicates(keep='first', inplace=True)

    ###Version 3, remove duplicates where the chr, start, and end for the query file are the same.
    if not boolKeepQueryTies:
        arUCEHeaders = arNewUCEHeaders[0:3]
        pdUCEsWithDistanceUnits.drop_duplicates(keep='first', inplace=True, subset=arUCEHeaders)

    return pdBreakpointsWithDistanceUnits, pdUCEsWithDistanceUnits


def reldist(pdAFile, pdBFile):
    btAFile = pandaToBedtool(pdAFile)
    btBFile = pandaToBedtool(pdBFile)
    dictReldist = btAFile.reldist(btBFile)
    #The columns are count, fraction, reldist, and total. To make a graph, you want to put the reldist on the x
    #axis, and the fraction on the y axis
    pdReldist = pd.DataFrame.from_dict(dictReldist)
    return pdReldist


def getSeriesWithBinsIndex(npHistDensityY, npHistDensityX):
    #Round the numbers
    npHistXRounded = np.around(npHistDensityX, decimals=3)

    #Get into list
    arHistXRounded = npHistXRounded.tolist()

    #Get rid of the last bin - there are always one more bin edges than there are numbers of bins, and you want
    #to make bin edges into lables (or indexes) for the counts, so you just forget the last one.
    #e.g. bin 0 to 0.1 is now just called bin 0, bin 0.49 to 0.5 is called bin 0.49
    arHistXRounded.pop()

    #Make into a series
    seriesHistDensity = pd.Series(npHistDensityY, index=arHistXRounded)

    return seriesHistDensity


def getNumpyHist(seriesDistances, intHistBins, tupXLim):
    #Get into list
    arDistances = seriesDistances.values.tolist()

    #Use np.histogram to split into 50 bins
    npHistCountsY, npHistCountsX = np.histogram(arDistances, bins=intHistBins, range=tupXLim, density=False)

    seriesCountsWithBins = getSeriesWithBinsIndex(npHistCountsY, npHistCountsX)

    return seriesCountsWithBins


def getIQR(pdRandomMerged):
    pdRandomMerged['upper_quartile'] = pdRandomMerged.apply(lambda row: row.quantile(q=0.75), axis=1)

    pdRandomMerged['median'] = pdRandomMerged.apply(lambda row: row.median(), axis=1)

    pdRandomMerged['lower_quartile'] = pdRandomMerged.apply(lambda row: row.quantile(q=0.25), axis=1)

    pdIQR = pdRandomMerged[['upper_quartile', 'median', 'lower_quartile']]

    return pdIQR


def doKS(seriesQueryRefReldist, seriesRandomReldistForKS):
    npQueryRefReldist = seriesQueryRefReldist.values
    npRandomReldist = seriesRandomReldistForKS.values
    ksStat, ksPval = stats.ks_2samp(npQueryRefReldist, npRandomReldist)
    print 'KS stat {0} KS P value {1}'.format(ksStat, ksPval)

    return ksStat, ksPval


def main():
    print 'Your query and ref files should have HEADERS.'
    print 'Your query and ref files should be ZERO based'
    print 'Your random files, if using, are expected to be ONE BASED!!!'

    args = get_args()
    strRefFilename = args.refFile
    strQueryFilename = args.queryFile

    strQueryLabel = args.queryLabel
    strRefLabel = args.refLabel

    intHistBins = args.bins
    tupXLim = args.xLims

    ###Version 3:
    boolKeepQueryTies = args.boolK

    #Get files as pandas
    pdRef = getDataIntoPandas(strRefFilename)
    pdQuery = getDataIntoPandas(strQueryFilename)

    #Get the distance from breakpoints to UCEs, and the distance from UCEs to breakpoints
    pdRefWithDistanceUnits, pdQueryWithDistanceUnits = mainGetDistances(pdRef, pdQuery, boolKeepQueryTies)

    #Save these results
    savePanda(pdRefWithDistanceUnits, 'Ref_{0}_with_distances_to_query_{1}.txt'.format(strRefLabel, strQueryLabel))
    savePanda(pdQueryWithDistanceUnits, 'Query_{0}_with_distances_to_ref_{1}.txt'.format(strQueryLabel, strRefLabel))

    #Get UCEs within 50kb
    intNoRefWithin50kb, pdRefWithin50kb = GetFeaturesWithinCertainDistance(pdRefWithDistanceUnits, 50)
    intNoRefWithin100kb, pdRefWithin100kb = GetFeaturesWithinCertainDistance(pdRefWithDistanceUnits, 100)

    savePanda(pdRefWithin50kb, 'Ref_{0}_with_query_{1}_regions_within_50kb.txt'.format(strRefLabel, strQueryLabel))
    savePanda(pdRefWithin100kb, 'Ref_{0}_with_query_{1}_regions_within_100kb.txt'.format(strRefLabel, strQueryLabel))

    intNoQueryWithin50kb, pdQueryWithin50kb = GetFeaturesWithinCertainDistance(pdQueryWithDistanceUnits, 50)
    intNoQueryWithin100kb, pdQueryWithin100kb = GetFeaturesWithinCertainDistance(pdQueryWithDistanceUnits, 100)

    savePanda(pdQueryWithin50kb, 'Query_{0}_with_ref_{1}_within_50kb.txt'.format(strQueryLabel, strRefLabel))
    savePanda(pdQueryWithin100kb, 'Query_{0}_with_ref_{1}_within_100kb.txt'.format(strQueryLabel, strRefLabel))

    #Get binned data for UCEs
    seriesQueryToRefKB = pdQueryWithDistanceUnits['Distance_Kb']
    ###Version 2: save to file
    savePanda(seriesQueryToRefKB, 'For_frequencies_query_{0}_to_ref_{1}_Kb_distances.txt'.format(strQueryLabel, strRefLabel))

    if args.filenames:
        arFilenames = get_Filenames_file(args.filenames)
        arPdRandoms = get_randoms_make_0based(arFilenames)

        arSeriesRandomBreakpointsDistBinned = []
        arSeriesRandomDistances = []

        for pdRandom in arPdRandoms:
            pdRandomsWithDistanceToRef, pdRefWithDistanceToRandoms = mainGetDistances(pdRandom, pdRef, boolKeepQueryTies)
            seriesDistancesKb = pdRandomsWithDistanceToRef['Distance_Kb']
            arSeriesRandomDistances.append(seriesDistancesKb)

            seriesRandomBreakpointsDistBinned = getNumpyHist(seriesDistancesKb, intHistBins, tupXLim)
            arSeriesRandomBreakpointsDistBinned.append(seriesRandomBreakpointsDistBinned)

        #Aggregate randoms
        pdRandomsDistanceCombined = pd.concat(arSeriesRandomBreakpointsDistBinned, axis=1)
        #Save this file
        pdRandomsDistanceCombined.to_csv('Random_to_ref_{0}_distance_as_counts.txt'.format(strRefLabel), sep="\t")

        #Get randoms IQR
        pdRandomsIQR = getIQR(pdRandomsDistanceCombined)


        #Get binned data for UCEs
        #seriesQueryToRefKB = pdQueryWithDistanceUnits['Distance_Kb']
        ###Version 2: save to file
        #savePanda(seriesQueryToRefKB, 'Query_{0}_to_ref_{1}_Kb_distances.txt'.format(strQueryLabel, strRefLabel))

        seriesQueryToRefDistBinned = getNumpyHist(seriesQueryToRefKB, intHistBins, tupXLim)

        #Put UCEs and randoms IQR together
        pdQueryAndRandomsBinnedForHist = pd.concat([seriesQueryToRefDistBinned, pdRandomsIQR], axis=1)

        pdQueryAndRandomsBinnedForHist.to_csv('Distances_query_{0}_to_ref_{1}_and_randomsIQR_to_ref_for_hist.txt'.format(strQueryLabel, strRefLabel), sep='\t')

        #Do KS
        seriesRandomDistances = pd.concat(arSeriesRandomDistances, ignore_index=True)
        ksStat, ksPval = doKS(seriesQueryToRefKB, seriesRandomDistances)
        with open('Ks_p_value', 'w') as fH:
            fH.write(str(ksPval))


        ##Version 2: for frequencies, save the series of distances in kb to separate files
        intNumberForSeries = 0
        arFilenamesForFrequencies = []
        for eachSeriesRandomDistances in arSeriesRandomDistances:
            eachFilename = 'For_frequences_random_to_ref_distances_{0}.txt'.format(intNumberForSeries)
            savePanda(eachSeriesRandomDistances, eachFilename)
            arFilenamesForFrequencies.append(eachFilename)
            intNumberForSeries += 1

        #Save these filenames to their own file to be read into script for binning frequencies
        with open('For_frequencies_filenames.txt', 'w') as fH:
            for strFile in arFilenamesForFrequencies:
                fH.write(strFile + '\n')

if __name__ == "__main__":
    main()