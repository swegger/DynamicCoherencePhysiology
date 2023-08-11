#!/usr/bin/env python

"""
SYNOPSIS

    findMaestroTrials SEARCHDIRECTORY LISTNAME TRIALNAME RESULTSFILE

DESCRIPTION

    Finds Maestro trials of the type TRIALNAME for each file in directory.

"""

import sys
import os
import csv

def main():
    queriedDir = sys.argv[1]
    listName = sys.argv[2]
    trialName = sys.argv[3]
    resultsFile = sys.argv[4]

    # Find directories
    dirList = listDirs(queriedDir)

    # create a list of trialTypes in each directory
    trialList(queriedDir,dirList,listName)

    # For each directory, determine number of maestro trials that match trialName
    trialNum = matchingTrials(queriedDir, dirList,listName,trialName)
    print(trialNum)
    # Export to files
    exportResults(queriedDir, dirList, trialNum, resultsFile)

    return 0


def listDirs(queriedDir):
    allfiles = sorted(os.listdir(queriedDir))
    out = []
    for item in allfiles:
        print(queriedDir + '/' + item)
        if os.path.isdir(queriedDir + '/' + item):
            out.append(item)

    return out

def trialList(queriedDir,directoryList,outputName):

    for filename in directoryList:
        f = os.path.join(queriedDir, filename)
        temp = f + '/' + outputName
        print(os.path.isfile(temp))
        if not os.path.isfile(temp):
            if os.path.isdir(f):
                a = 'find '
                b = f
                c = "/* -exec head -1 \{\} \; > "
                d = f
                e = outputName
                bashCommand = a + b + c + d + '/' + e
    #            print(bashCommand)
                os.system(bashCommand)

    return 0

def matchingTrials(queriedDir, directoryList, outputName, trialName):
    trialName_ = trialName + '_'
    trialNum = []
    for fileIdx, filename in enumerate(directoryList):
        tempFile = open(queriedDir + '/' + filename + '/' + outputName, 'r', errors='ignore')
        Lines = tempFile.readlines()

        count = 0
        for lines in Lines:
            print(lines)
            if trialName_ in lines:
                count += 1

        print(count)
        trialNum.append(count)

    return trialNum

def exportResults(queriedDir, directoryList, trialNumber, outputFile):

      with open (queriedDir + '/' + outputFile + '.csv', "w", newline = '') as csvfile:
          header = ['File', 'TrialN']
          my_writer = csv.DictWriter(csvfile, fieldnames = header)
          my_writer.writeheader()
          for idx, temp in enumerate(trialNumber):
              my_writer.writerow({'File' : directoryList[idx],
                                  'TrialN' : trialNumber[idx]})

if __name__ == '__main__':
    raise SystemExit(main())
