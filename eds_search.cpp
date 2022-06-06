//#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "globals.hpp"
#include "functions.hpp"
#include "bndm.hpp"
#include "translator.hpp"

#include "cedsm.hpp"
#include "helpers.hpp"
#include "params.hpp"
#include "parsing.hpp"
#include "sopang.hpp"

using namespace sopang;
using namespace std;

struct SegmentData
{
    const string *const *segments;
    int nSegments; // Number of segments.
    const int *segmentSizes; // Size of each segment (number of variants).
};

string readInputText(const string fName, double &textSizeMB)
{
    string text = helpers::readFile(fName);
    cout << "Read file: " << fName << endl;

    textSizeMB = text.size() / 1'000'000.0;

    /* if (params.dumpToFile)
    {
        string outStr = params.inTextFile + " " + to_string(textSizeMB) + " ";
        helpers::dumpToFile(outStr, params.outFile, false);
    } */

    cout << boost::format("Read input text, #chars = %1%, MB = %2%") % text.size() % textSizeMB << endl;
    return text;
}

vector<vector<Sopang::SourceSet>> readSources(int nSegments, const int *segmentSizes, int &sourceCount, string sourceFilePath)
{
    string sourcesStr = helpers::readFile(sourceFilePath);
    cout << "Read file: " << sourceFilePath << endl;

    vector<vector<Sopang::SourceSet>> sources;
    cout << "Parsing sources..." << endl;
    sources = parsing::parseSources(sourcesStr, sourceCount);

    //cout << boost::format("Parsed sources for non-deterministic #segments = %1%, #chars = %2%")
    //    % sources.size() % sourcesStr.size() << endl;

    cout << "Source count = " << sourceCount << endl;

    if (sources.empty())
    {
        throw runtime_error("cannot run for empty sources");
    }

    // We check whether the source counts match the segments in text.
    size_t sourceIdx = 0;

    for (int segmentIdx = 0; segmentIdx < nSegments; ++segmentIdx)
    {
        if (segmentSizes[segmentIdx] == 1)
        {
            continue;
        }

        if (static_cast<size_t>(segmentSizes[segmentIdx]) != sources[sourceIdx].size())
        {
            throw runtime_error("source segment variant count does not match text segment variant count, source segment index = "
                + to_string(sourceIdx));
        }

        // We check whether all sources are present for the current segment.
        Sopang::SourceSet sourcesForSegment(sourceCount);

        for (const Sopang::SourceSet &sourcesForVariant : sources[sourceIdx])
        {
            sourcesForSegment |= sourcesForVariant;
        }

        if (sourcesForSegment.count() != sourceCount)
        {
            throw runtime_error("not all sources are present for segment: " + to_string(sourceIdx));
        }

        sourceIdx += 1;

        if (sourceIdx > sources.size())
        {
            throw runtime_error("there are fewer source segments than non-deterministic segments in text");
        }
    }

    if (sourceIdx < sources.size())
    {
        throw runtime_error("there are more source segments than non-deterministic segments in text");
    }

    cout << "Sanity check for sources passed" << endl;
    return sources;
}

std::string readFile(const std::string &filePath)
{
    using namespace std;

    ifstream inStream(filePath);

    if (!inStream)
    {
        throw runtime_error("failed to read file (insufficient permisions?): " + filePath);
    }

    return static_cast<stringstream const&>(stringstream() << inStream.rdbuf()).str();
}

int main(int argc, char *argv[]) {

    int LOOPS = atoi(argv[2]);
    size_t pattLen = (unsigned int) atoi(argv[3]), patt2Len = 0;
    std::string outFile = "";
    const char *algorithm = argv[4];
    unsigned char patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH + 1];
    memset(patterns, 0, MAX_DNA_PATTERNS*(MAX_PATTERN_LENGTH + 1));

    if (argc >= 6 && *argv[5]) {
        pattLen = strnlen(argv[5], MAX_PATTERN_LENGTH + 1);
        memcpy(patterns[0], argv[5], pattLen);
    }
    /* if (argc >= 7 && *argv[6]) {
        patt2Len = strnlen(argv[6], MAX_PATTERN_LENGTH + 1);
        strncpy(patterns[1], argv[6], pattLen);
    } */
    if (pattLen && patt2Len && pattLen!=patt2Len) {
        fprintf(stderr, "Two patterns supplied, but length is not matching!\n");
        return EXIT_FAILURE;
    }
    if (argc >= 8 && *argv[7]) {
        outFile = argv[7];
    }

    int eds_size, teds_size;
    unsigned char *eds = readInputFile(argv[1], &eds_size);
    unsigned char *teds = translate(eds, eds_size, &teds_size);
    free(eds);

//    writeOutputFile("/tmp/eds_search.teds", &teds, teds_size);

    //init_IUPAC_SYMBOLS_TO_BASES();
    //init_AA_TO_COMPR_IUPAC_SYMBOLS();

    srand(time(NULL));

    if (strcasecmp(algorithm, "bndm") == 0) {
        bndm_eds_run(teds, teds_size, patterns[0], pattLen, LOOPS);
    }
    else if (strcasecmp(algorithm, "cedsm") == 0) {
        // parse text
        double textSizeMB;
        const string text = readInputText(argv[1], textSizeMB);
        cout << "Parsing segments..." << endl;

        int nSegments = 0;
        int *segmentSizes = nullptr;

        const string *const *segments = parsing::parseTextArray(text, &nSegments, &segmentSizes);
        cout << "Parsed #segments = " << nSegments << endl;
        SegmentData segmentData{ segments, nSegments, segmentSizes };

        // parse sources
        int sourceCount = 0;
        Sopang::SourceMap sourceMap;
        const vector<vector<Sopang::SourceSet>> sources = readSources(nSegments, segmentSizes, sourceCount, argv[6]);
        sourceMap = parsing::sourcesToSourceMap(nSegments, segmentSizes, sources);

        cedsm_eds_run(  segmentData.segments,
                        segmentData.nSegments,
                        segmentData.segmentSizes,
                        patterns[0], pattLen, sourceMap, sourceCount, LOOPS, textSizeMB, outFile);
    } else {
        fprintf(stderr, "Invalid algorithm option!\n");
        free(teds);
        return EXIT_FAILURE;
    }
    free(teds);
    return EXIT_SUCCESS;
}