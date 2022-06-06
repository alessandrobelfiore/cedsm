"""
Checks whether the result (number of matches) from sopang is the same as the desired values.


Usage: count_results.py <file1> <file2> [options]

Arguments:
  <file1>       path to the first file (ceds)
  <file2>       path to the second file (sopang)
"""

import sys
from docopt import docopt

def readMatchesFromFile(inFile):
    with open(inFile, "r") as f:
        data = f.read()

    lines = data.split("\n")
    lines = [l for l in lines if l]

    #matches = []
    matches = 0

    for l in lines:
        if l.startswith("#results"):
            #matches += [int(l.split()[-1])]
            matches += int(l.split()[-1])

    print("Got #matches = {0} from file".format(matches))
    return matches

""" def readMatchesFromCmdLine():
    matches = sys.argv[1].split()
    matches = [int(n) for n in matches]

    print("Got #matches = {0} from cmd-line".format(len(matches)))
    return matches """

def main():
    args = docopt(__doc__, version="0.1.0")
    pInFile = args["<file1>"]
    pInFile2 = args["<file2>"]

    fileMatches = readMatchesFromFile(pInFile)
    fileMatches2 = readMatchesFromFile(pInFile2)
    #cmdMatches = readMatchesFromCmdLine()
    difference = fileMatches2 - fileMatches
    percentage = (difference * 100) / fileMatches2

    #print("OK (matches = {0})".format(fileMatches))
    #print("OK (matches = {0})".format(fileMatches2))
    print("Difference: {0}".format(difference))
    print("Percentage: {0}".format(percentage))
    #else:
    #    print("ERROR: mismatches")

if __name__ == "__main__":
    main()
