"""
Script to implement Anderson darling test for comparison of distributions. Takes 2 files that are generated by
distances2.py script: Random_to_ref_distances_as_list.txt and a file with name in form
Query_UCEs_with_distances_to_ref_Ear_eye_nervous_Pooledbreakpoints_1kbcluster.txt

Ruth McCole
April 3rd 2017

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
from scipy import stats


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument('-q', "--queryFile", type=str,
                        help='Filename for query to ref distances')
    parser.add_argument('-r', '--randomFile', type=str,
                        help='Filename for random to ref distances')
    parser.add_argument('-l', '--limit', type=float,
                        help='Limit for lengths to cutoff for distribution comparison')
    parser.add_argument('-o', '--outputStem', type=str, default='Test',
                        help='Output stem')

    return parser.parse_args()


def getDataIntoPandas(strFilename):
    #NB: this assumes your files DO have headers
    print 'Getting data from {0}'.format(strFilename)
    pdData = pd.read_csv(strFilename, sep='\t')
    return pdData


def getDataIntoPandasNoHeader(strFilename, arHeaders):
    #This is for files without headers
    pdData = pd.read_csv(strFilename, sep='\t', header=None, names=arHeaders)
    return pdData

def savePanda(pdData, strFilename):
    pdData.to_csv(strFilename, sep='\t', index=False)


def applyLimit(pdDistances, floatLimit, strColumnName):
    #Now have to make this work as a series, or make randoms a 1 column dataframe
    BoolPdLimit = pdDistances[strColumnName] < floatLimit
    pdDistancesLimit = pdDistances[BoolPdLimit]
    print 'Limit of {0}kb applied'.format(floatLimit)

    return pdDistancesLimit


def getAD(pdQuery, pdRef):
    arQuery = pdQuery.values.flatten().tolist()
    arRef = pdRef.values.flatten().tolist()

    tupAD = stats.anderson_ksamp([arQuery, arRef])

    return tupAD


def processTupAD(tupAD, strOutputStem):
    print 'AD test stat is {0}'.format(str(tupAD[0]))
    print 'AD critical values at 25%, 10%, 5%, 2.5%, 1% respecitively are:'
    print tupAD[1]
    arCritical = tupAD[1]
    print 'Approx p-value is {0}'.format(str(tupAD[2]))

    with open('AD_result_{0}.txt'.format(strOutputStem), 'a') as fH:
        fH.write('AD test stat is {0}'.format(str(tupAD[0])) + "\n")
        fH.write('AD critical values at 25%, 10%, 5%, 2.5%, 1% respectively are:' +'\n')
        fH.write(str(tupAD[1]) + '\n')
        fH.write('Approx p-value is {0}'.format(str(tupAD[2])))


def main():
    args = get_args()

    strOutputStem = args.outputStem
    pdQueryFull = getDataIntoPandas(args.queryFile)

    pdQueryDistances = pdQueryFull[['Distance_Kb']]

    seriesRandomDistances = pd.read_csv(args.randomFile, squeeze=True)
    pdRandomDistances = getDataIntoPandasNoHeader(args.randomFile, ['Random_distances_Kb'])


    if args.limit:
        floatLimit = args.limit
        print 'limit of {0}kb'.format(str(floatLimit))
        pdQueryLimit = applyLimit(pdQueryDistances, floatLimit, 'Distance_Kb')
        pdRandomsLimit = applyLimit(pdRandomDistances, floatLimit, 'Random_distances_Kb')

        #Save for future ref
        savePanda(pdQueryLimit, 'Query_limit_to_{0}kb.txt'.format(str(floatLimit)))
        savePanda(pdRandomsLimit, 'Randoms_limited_to_{0}kb.txt'.format(str(floatLimit)))

        tupAD = getAD(pdQueryLimit, pdRandomsLimit)
        processTupAD(tupAD, strOutputStem)

    else:
        tupAD = getAD(pdQueryDistances, pdRandomDistances)
        processTupAD(tupAD, strOutputStem)


if __name__ == "__main__":
    main()




