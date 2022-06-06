#include "sopang.hpp"

#include <algorithm>
#include <cassert>

#include <iostream>
#include <cmath>
#include <list>
#include <cstring>

using namespace std;
using namespace sopang;

unordered_set<int> Sopang::matchCEDS(const string *const *segments,
    int nSegments,
    const int *segmentSizes,
    const Sopang::SourceMap &sourceMap,
    int sourceCount,
    const string &pattern)
{
    assert(nSegments > 0 && pattern.size() > 0 && pattern.size() <= wordSize);

    fillPatternMaskBuffer(pattern);

    const uint64_t hitMask = (0x1ULL << (pattern.size() - 1));
    uint64_t D = allOnes;

    unordered_set<int> res;
    // Search phase
    for (int iS = 0; iS < nSegments; ++iS)  // for every segment
    {
        assert(segmentSizes[iS] > 0 && static_cast<size_t>(segmentSizes[iS]) <= dBufferSize);

        for (int iD = 0; iD < segmentSizes[iS]; ++iD)   // for every variant present in the segment
        {
            dBuffer[iD] = D;
            
            for (int i = 0; i < pattern.size(); ++i)
            {
                if ((dBuffer[iD] & indexHitMask[i]) == 0x0ULL)
                    *(sBuffer[iD][i]) = *(S[i]);
            }

            if (segments[iS][iD].size() == 0)   // if the variant is empty, we still need to update the sources
            {
                if (sourceMap.find(iS) != sourceMap.end())  // if it's not a solid segment
                {
                    for (size_t i = 0; i < pattern.size(); ++i) // for every char in pattern
                    {
                        if ((dBuffer[iD] & indexHitMask[i]) == 0x0ULL)
                        {
                            *(sBuffer[iD][i]) = *(sBuffer[iD][i]) & sourceMap.at(iS)[iD];
                        }
                    }
                }
            }

            for (size_t iC = 0; iC < segments[iS][iD].size(); ++iC) // for every char in the string
            {
                const char c = segments[iS][iD][iC];

                assert(c > 0 && static_cast<unsigned char>(c) < maskBufferSize);
                assert(alphabet.find(c) != string::npos);

                dBuffer[iD] <<= 1;
                dBuffer[iD] |= maskBuffer[static_cast<unsigned char>(c)];
                // shift sources
                shiftSources(sBuffer[iD], pattern.size());
                // we combine the available sources on the active bits
                if (sourceMap.find(iS) != sourceMap.end())
                {
                    for (size_t i = 0; i < pattern.size(); ++i) 
                    {
                        if ((dBuffer[iD] & indexHitMask[i]) == 0x0ULL)
                        {
                                #ifdef DEBUG
                                std::cout << "ANDing: " << 
                                sBuffer[iD][i].buffer[0] << " with " << sourceMap.at(iS)[iD].buffer[0] << std::endl;
                                #endif
                                *(sBuffer[iD][i]) = *(sBuffer[iD][i]) & sourceMap.at(iS)[iD];
                                #ifdef DEBUG
                                std::cout << "Resulting in: " << sBuffer[iD][i].buffer[0] << " index is: " <<
                                iS << " " << iD << " " << iC << std::endl;
                                #endif
                        }
                    }
                }
                // Match occurred. Note: we still continue in order to fill the whole d-buffer.
                if ((dBuffer[iD] & hitMask) == 0x0ULL)
                {
                    #ifdef DEBUG
                    std::cout << "found a match! " << " index is: " <<
                            iS << " " << iD << " " << iC << std::endl;
                    for (size_t i = 0; i < sBuffer[iD][0].bufferSize; ++i)
                    {
                        std::cout << sBuffer[iD][0].buffer[i] << std::endl;
                    }
                    #endif

                    if (sBuffer[iD][pattern.size() - 1]->any()) // if at least one source is present
                    {
                        res.insert(iS);
                    }
                }
            }
        }
        // End of segment, we need to merge results from variants
        D = dBuffer[0];
        for (size_t i = 0; i < pattern.size(); ++i) // for every char
        {
            if ((dBuffer[0] & indexHitMask[i]) == 0x0ULL) // if the bit is active (== 0)
            {
                *(S[i]) = *(sBuffer[0][i]);
            } 
            else 
            {
                S[i]->reset();
            }
        }

        for (int iD = 1; iD < segmentSizes[iS]; ++iD)   // for every variation
        {
            // As a join operation we want to preserve 0s (active states):
            // a match can occur in any segment alternative.
            D &= dBuffer[iD];

            // for every active bit we add all possible sources from the active bits
            for (size_t i = 0; i < pattern.size(); ++i) // for every char
            { 
                if ((dBuffer[iD] & indexHitMask[i]) == 0x0ULL) // if the bit is active (== 0)
                {
                    #ifdef DEBUG
                    std::cout << "ORing: " << 
                        sBuffer[iD][i].buffer[0] << 
                        " with " << S[i].buffer[0] << std::endl;
                    #endif
                    *(S[i]) |= *(sBuffer[iD][i]);
                    #ifdef DEBUG
                    std::cout << "Resulting in: " << S[i].buffer[0] << " index is: " <<
                        iS << " " << iD << " " << i << std::endl;
                    #endif
                }
            }
        }
    }
    return res;
}