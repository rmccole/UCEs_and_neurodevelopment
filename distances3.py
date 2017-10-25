"""
This script takes outputs from distances2.py and bins the data, producing a final results table
ready to be plotted as a histogram.

Ruth McCole
February 6th 2017

Copyright 2017 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import argparse
import pandas as pd
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument('-q', "--queryToRef", type=str,
                        help='Filename for the query to ref distances')
    parser.add_argument('-f', '--filenames', type=str,
                        help='File containing filenames of files containing query to random distances')
    parser.add_argument('-rl', '--refLabel', type=str,
                        help='Label for reference')
    parser.add_argument('-ql', '--queryLabel', type=str,
                        help='Label for query')
    parser.add_argument('-b', '--bins', type=int, default=50,
                        help='Number of bins to aggregate into')
    parser.add_argument('-x', '--xLims', type=int, nargs=2, default=(0, 20000),
                        help='Upper and lower limits of bins for data aggregation, supply two integers, separated by'
                             'a space')

    return parser.parse_args()


def get_Filenames_file(strFilenames):
    """
    Obtain a list of the filenames of the files to process
    """
    with open(strFilenames, 'r') as Filenames_filehandle:
        aFiles = [line.strip() for line in Filenames_filehandle]
    return aFiles


def getDataIntoPandasNoHeader(strFilename, arHeaders):
    #This is for files without headers
    #print 'Getting data from {0}'.format(strFilename)
    pdData = pd.read_csv(strFilename, sep='\t', header=None, names=arHeaders)
    return pdData


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


def getIQRAndRange(pdRandomMerged):
    pdRandomMerged['upper_quartile'] = pdRandomMerged.apply(lambda row: row.quantile(q=0.75), axis=1)

    pdRandomMerged['median'] = pdRandomMerged.apply(lambda row: row.median(), axis=1)

    pdRandomMerged['lower_quartile'] = pdRandomMerged.apply(lambda row: row.quantile(q=0.25), axis=1)

    pdRandomMerged['min'] = pdRandomMerged.apply(lambda row: row.min(), axis=1)

    pdRandomMerged['max'] = pdRandomMerged.apply(lambda row: row.max(), axis=1)

    pdIQRAndRange = pdRandomMerged[['upper_quartile', 'median', 'lower_quartile', 'min', 'max']]

    return pdIQRAndRange


def main():
    args = get_args()

    intHistBins = args.bins
    tupXLim = args.xLims
    intXMax = tupXLim[-1]

    strQueryLabel = args.queryLabel
    strRefLabel = args.refLabel

    #Get query to ref distances
    pdQueryToRefDist = getDataIntoPandasNoHeader(args.queryToRef, ['Query_to_ref'])
    seriesQueryToRef = pdQueryToRefDist['Query_to_ref']

    #Bin query to ref distances
    seriesQueryCountsWithBins = getNumpyHist(seriesQueryToRef, intHistBins, tupXLim)
    seriesQueryCountsWithBins.rename('{0}_to_{1}_distance_counts'.format(strQueryLabel, strRefLabel), inplace=True)

    #Get list of series for randoms
    arRandomFilenames = get_Filenames_file(args.filenames)

    arSeriesRandomBinned = []
    for strRandomFilename in arRandomFilenames:
        pdRandomToRefDist = getDataIntoPandasNoHeader(strRandomFilename, ['Random_to_ref'])
        seriesRandomToRef = pdRandomToRefDist['Random_to_ref']
        seriesRandomCountsWithBins = getNumpyHist(seriesRandomToRef, intHistBins, tupXLim)
        arSeriesRandomBinned.append(seriesRandomCountsWithBins)

    #Aggregate randoms, forms a dataframe 1000 columns wide
    pdRandomsDistanceCombined = pd.concat(arSeriesRandomBinned, axis=1)
    #Save this file
    pdRandomsDistanceCombined.to_csv('Random_to_ref_{0}_distance_as_counts_{1}_bins_{2}_xlim.txt'.format(strRefLabel, intHistBins, intXMax), sep="\t")

    #Get randoms IQR and range
    pdRandomsIQRAndRange = getIQRAndRange(pdRandomsDistanceCombined)

    #Put UCEs and randoms IQR together
    pdQueryAndRandomsBinnedForHist = pd.concat([seriesQueryCountsWithBins, pdRandomsIQRAndRange], axis=1)

    pdQueryAndRandomsBinnedForHist.to_csv('Distances_{0}_to_{1}_and_randomsIQRRange_to_ref_for_hist_for_{2}_bins_{3}_xlim.txt'.format(strQueryLabel, strRefLabel, intHistBins, intXMax), sep='\t')


if __name__ == "__main__":
    main()
