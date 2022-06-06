#include "sopang.hpp"

#include <algorithm>
#include <cassert>

#include <iostream>
#include <cmath>
#include <list>
#include <cstring>
/* #define DEBUG */

#define SIGMA 256       // Max size of the alphabet

using namespace std;

namespace sopang
{

Sopang::Sopang(const std::string &alphabet, int sourceCount)
    :alphabet(alphabet)
{
    dBuffer = new uint64_t[dBufferSize];
    initCounterPositionMasks();
    initIndexHitMask();
    initSourcesVector(sourceCount);
    Sopang::sourceCount = sourceCount;
}

void Sopang::initIndexHitMask()
{
    for (size_t i = 0; i < wordSize; ++i)
    {
        Sopang::indexHitMask[i] = 0x1ULL << i;
    }
}

void Sopang::initSourcesVector(int sourceCount)
{
    for (size_t i = 0; i < wordSize; ++i)
    {
        S[i] = new Sopang::SourceSet (sourceCount);
        S[i]->set();
    }
    for (int i = 0; i < sourceCount; ++i) // for every variant
    {
        for (size_t j = 0; j < wordSize; ++j) // for every character
        {
            sBuffer[i][j] = new Sopang::SourceSet (sourceCount);
            sBuffer[i][j]->set();
        }
    }
}

Sopang::~Sopang()
{
    delete[] dBuffer;
    for (size_t i = 0; i < wordSize; ++i) 
    {
        delete(S[i]);
    }
    
    for (int i = 0; i < sourceCount; ++i) // for every variant
    {
        for (size_t j = 0; j < wordSize; ++j) // for every character
        {
            delete(sBuffer[i][j]);
        }
    }
}

void Sopang::shiftSources(Sopang::SourceSet* sources [], int size)   // rotates 1 position right
{
    Sopang::SourceSet* tmp = sources[size - 1];
    for (int i = size - 1; i > 0; --i)
    {
        sources[i] = sources[i - 1];
    }
    sources[0] = tmp;
    sources[0]->set();
}

std::unordered_set<int> Sopang::matchBNDM_CEDS(const std::string *const *segments,
		int nSegments,
        const int *segmentSizes,
        const Sopang::SourceMap &sourceMap,
        int sourceCount,
        const std::string &pattern) 
{
	uint64_t m = pattern.size();
	unsigned int B[SIGMA];
	unsigned int maskBuffer[SIGMA];
	unsigned int* R = (unsigned int*)malloc(m*sizeof(unsigned int)); //Stores starting register values for change from BNDM to SA. Basically one-set-bit on different positions 0 to m-1.
	int i, j, hitMask, D, D2, last, R1, R2;
	char c;
	bool isSolid;

	// Sources Preprocessing
	Sopang::SourceSet* sLast[wordSize];

	for (size_t i = 0; i < wordSize; ++i)
	{
		sLast[i] = new Sopang::SourceSet (sourceCount);
	}

	// BNDM Preprocessing
	for (i = 0; i < SIGMA; i++) B[i] = 0;
	hitMask = 1;
	
	for (i = m - 1; i >= 0; i--) {
		B[pattern[i]] |= hitMask;
		R[i] = hitMask;
		hitMask <<= 1;
	}
	hitMask >>= 1;
	// SA Preprocessing
	for (i = 0; i < SIGMA; ++i) 
		maskBuffer[i] = 0;
	for (i = 0, j = 1; i < m; ++i, j <<= 1)
		maskBuffer[pattern[i]] |= j;
	
	unsigned int elementLength;
	unsigned int elementStart, elementEnd;

	R1 = R2 = D2 = 0;
	uint64_t* dBuffer = new uint64_t[sourceCount];
	std::unordered_set<int> res;

	// Searching phase
	for (int iS = 0; iS < nSegments; ++iS)  // for every segment
    {
		R1 = R2;
		R2 = 0;
		isSolid = sourceMap.find(iS) != sourceMap.end() ? false : true;
		for (int iD = 0; iD < segmentSizes[iS]; ++iD)   // for every variant present in the segment
		{
			elementEnd = elementLength = segments[iS][iD].size();
			j = 0;
			// Re-initialize dBuffer
			dBuffer[iD] = R1;
			// Re-initialize sBuffer
			for (int i = 0; i < m; ++i)
            {
                if (dBuffer[iD] & indexHitMask[i])
                    *(sBuffer[iD][i]) = *(S[i]);
            }
			// Check if the variant is empty
			if (segments[iS][iD].size() == 0)   // if the variant is empty, we still need to update the sources
            {
                if (!isSolid)  // if it's not a solid segment
                {
                    for (size_t i = 0; i < m; ++i) // for every char in pattern
                    {
                        if (dBuffer[iD] & indexHitMask[i])
                        {
                            *(sBuffer[iD][i]) = *(sBuffer[iD][i]) & sourceMap.at(iS)[iD];
                        }
                    }
                }
            }
			// Process first m bases with SA
			while (j < m && j < elementEnd)
			{
				c = segments[iS][iD][j];
				dBuffer[iD] = ((dBuffer[iD] << 1) | 1) & maskBuffer[static_cast<unsigned char>(c)];
				shiftSources(sBuffer[iD], m);
				// Intersect sources
				if (!isSolid)
                {
                    for (size_t i = 0; i < m; ++i)
                    {
                        if (dBuffer[iD] & indexHitMask[i])
                        {
                            *(sBuffer[iD][i]) = *(sBuffer[iD][i]) & sourceMap.at(iS)[iD];
                        }
                    }
                }
				// Check match
				if (dBuffer[iD] & hitMask)
				{
					if (sBuffer[iD][m - 1]->any()) // if at least one source is present
                    {
						res.insert(iS);
                    }
				}
				j++;
			}
			if (elementLength < m)
				continue;

			// Perform BNDM search.
			j = 0;
			while (j + m - 1 < elementEnd)
			{
				D2 = 0;
				last = m;
				i = m - 1;
				dBuffer[iD] = ~0;
				// Re-initialize sources
				for (int k = 0; k < m; ++k)
				{
					sBuffer[iD][k]->set();
				}
				while (i >= 0 && dBuffer[iD] != 0) 
				{
					c = segments[iS][iD][j + i];
					dBuffer[iD] &= B[static_cast<unsigned char>(c)];
					// Intersect sources
					if (!isSolid)
					{
						for (unsigned int i = 0; i < m; ++i)
						{
							if (dBuffer[iD] & indexHitMask[m - i - 1])
							{
								*(sBuffer[iD][i]) = *(sBuffer[iD][i]) & sourceMap.at(iS)[iD];
							}
						}
					}
					i--;
					// Found a prefix
					if (dBuffer[iD] != 0  && (dBuffer[iD] & hitMask) != 0)
					{
						if (i >= 0)
						{
							last = i + 1;
							D2 |= R[last];
							*(sLast[m - (last + 1)]) = *(sBuffer[iD][m - (last + 1)]);
						}
						else
						{
							// Check match
							if (sBuffer[iD][m - 1]->any()) // if at least one source is present
							{
								res.insert(iS);
							}
						}
					}
					dBuffer[iD] <<= 1;
				}
				j += last;
			}
			//Perform SA search at the end of the element.
			j += m - last; //Moving j pointer to the initial position for SA.
			dBuffer[iD] = D2; //Setting the initial value to SA register.
			for (int i = 0; i < m; ++i)
			{
				if (dBuffer[iD] & indexHitMask[i])
				{
					*(sBuffer[iD][i]) = *(sLast[i]);
				}
				else 
					sBuffer[iD][i]->set();
			}
			while (j < elementEnd)
			{
				c = segments[iS][iD][j];
				dBuffer[iD] = ((dBuffer[iD] << 1) | 1) & maskBuffer[static_cast<unsigned char>(c)];
				
				shiftSources(sBuffer[iD], m);
				// Intersect sources
				if (sourceMap.find(iS) != sourceMap.end())
				{
					for (unsigned int i = 0; i < m; ++i) 
					{
						if (dBuffer[iD] & indexHitMask[i])
						{
							*(sBuffer[iD][i]) = *(sBuffer[iD][i]) & sourceMap.at(iS)[iD];
						}
					}
				}
				// Checking match
				if (dBuffer[iD] & hitMask)
				{
					if (sBuffer[iD][m - 1]->any()) // if at least one source is present
                    {
						res.insert(iS);
                    }
				}
				j++;
			}
		}
		// Merging intermediate results 
		R2 = dBuffer[0];
		for (unsigned int i = 0; i < m; ++i) // for every char
        {
            if (dBuffer[0] & indexHitMask[i]) // if the bit is active (== 1)
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
            R2 |= dBuffer[iD];
            // for every active bit we add all possible sources from the active bits
            for (unsigned int i = 0; i < m; ++i) // for every char
            {
                if (dBuffer[iD] & indexHitMask[i]) // if the bit is active (== 0)
                {
                    *(S[i]) |= *(sBuffer[iD][i]);
                }
            }
        }
	}
	return res;
}

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

unordered_set<int> Sopang::match(const string *const *segments,
    int nSegments,
    const int *segmentSizes,
    const string &pattern)
{
    assert(nSegments > 0 and pattern.size() > 0 and pattern.size() <= wordSize);

    fillPatternMaskBuffer(pattern);

    const uint64_t hitMask = (0x1ULL << (pattern.size() - 1));
    uint64_t D = allOnes;

    unordered_set<int> res;

    for (int iS = 0; iS < nSegments; ++iS)
    {
        assert(segmentSizes[iS] > 0 and static_cast<size_t>(segmentSizes[iS]) <= dBufferSize);

        for (int iD = 0; iD < segmentSizes[iS]; ++iD)
        {
            dBuffer[iD] = D;

            for (size_t iC = 0; iC < segments[iS][iD].size(); ++iC)
            {
                const char c = segments[iS][iD][iC];

                assert(c > 0 and static_cast<unsigned char>(c) < maskBufferSize);
                assert(alphabet.find(c) != string::npos);

                dBuffer[iD] <<= 1;
                dBuffer[iD] |= maskBuffer[static_cast<unsigned char>(c)];

                // Match occurred. Note: we still continue in order to fill the whole d-buffer.
                if ((dBuffer[iD] & hitMask) == 0x0ULL)
                {
                    res.insert(iS);
                    #ifdef DEBUG
                    std::cout << "inserting result at: " << iS << " " 
                    << iD << " " << iC << std::endl;
                    #endif
                }
            }
        }

        D = dBuffer[0];

        for (int iD = 1; iD < segmentSizes[iS]; ++iD)
        {
            // As a join operation we want to preserve 0s (active states):
            // a match can occur in any segment alternative.
            D &= dBuffer[iD];
        }
    }

    return res;
}

unordered_set<int> Sopang::matchApprox(const string *const *segments,
    int nSegments,
    const int *segmentSizes,
    const string &pattern,
    int k)
{
    assert(nSegments > 0 and pattern.size() > 0 and pattern.size() <= maxPatternApproxSize);
    assert(k > 0);

    unordered_set<int> res;

    fillPatternMaskBufferApprox(pattern);

    // This is the initial position of each counter, after k + 1 errors the most significant bit will be set.
    const uint64_t counterMask = 0xFULL - k;
    // Hit mask indicates whether the most significant bit in the last counter is set.
    const uint64_t hitMask = (0x1ULL << ((pattern.size() * saCounterSize) - 1));

    uint64_t D = 0x0ULL;

    for (size_t i = 0; i < pattern.size(); ++i)
    {
        // We initialize the state with full counters which allow us to effectively start matching
        // after m characters.
        D |= (saFullCounter << (i * saCounterSize));
    }

    for (int iS = 0; iS < nSegments; ++iS)
    {
        assert(segmentSizes[iS] > 0 and static_cast<size_t>(segmentSizes[iS]) <= dBufferSize);

        for (int iD = 0; iD < segmentSizes[iS]; ++iD)
        {
            dBuffer[iD] = D;

            for (size_t iC = 0; iC < segments[iS][iD].size(); ++iC)
            {
                const char c = segments[iS][iD][iC];

                assert(c > 0 and static_cast<unsigned char>(c) < maskBufferSize);
                assert(alphabet.find(c) != string::npos);

                dBuffer[iD] <<= saCounterSize;
                dBuffer[iD] += counterMask;

                dBuffer[iD] += maskBuffer[static_cast<unsigned char>(c)];

                if ((dBuffer[iD] & hitMask) == 0x0ULL)
                {
                    res.insert(iS);
                }
            }
        }

        D = 0x0ULL;

        // As a join operation, we take the minimum (the most promising alternative) from each counter.
        for (size_t i = 0; i < pattern.size(); ++i)
        {
            uint64_t min = (dBuffer[0] & counterPosMasks[i]);

            for (int iD = 1; iD < segmentSizes[iS]; ++iD)
            {
                uint64_t cur = (dBuffer[iD] & counterPosMasks[i]);
                
                if (cur < min)
                {
                    min = cur;
                }
            }

            D |= min;
        }
    }

    return res;
}

namespace
{

bool verifyMatch(const string *const *segments,
    const int *segmentSizes,
    const Sopang::SourceMap &sourceMap,
    int sourceCount,
    const string &pattern,
    int matchIdx,
    const pair<int, int> &match)
{
    using SourceSet = Sopang::SourceSet;
    const int patternCharIdx = static_cast<int>(pattern.size()) - match.second - 2;

    if (patternCharIdx < 0) // The match is fully contained within a single segment.
        return true;

    vector<pair<SourceSet, int>> leaves;
    SourceSet rootSources(sourceCount);

    if (sourceMap.count(matchIdx) > 0)
    {
        rootSources = sourceMap.at(matchIdx)[match.first];
    }
    else
    {
        rootSources.set();
    }

    leaves.emplace_back(move(rootSources), patternCharIdx);
    int segmentIdx = static_cast<int>(matchIdx - 1);

    while (not leaves.empty() and segmentIdx >= 0)
    {
        if (segmentSizes[segmentIdx] == 1)
        {
            for (auto &leaf : leaves)
            {
                leaf.second -= segments[segmentIdx][0].size();

                if (leaf.second < 0)
                    return true;
            }
        }
        else
        {
            assert(sourceMap.count(segmentIdx) > 0 and sourceMap.at(segmentIdx).size() == static_cast<size_t>(segmentSizes[segmentIdx]));

            vector<pair<SourceSet, int>> newLeaves;
            newLeaves.reserve(leaves.size() * segmentSizes[segmentIdx]);

            for (const auto &leaf : leaves)
            {
                for (int variantIdx = 0; variantIdx < segmentSizes[segmentIdx]; ++variantIdx)
                {
                    const SourceSet &variantSources = sourceMap.at(segmentIdx)[variantIdx];

                    if (segments[segmentIdx][variantIdx].empty())
                    {
                        const SourceSet newSources = (variantSources & leaf.first);

                        if (newSources.any())
                        {
                            newLeaves.emplace_back(move(newSources), leaf.second);
                        }
                    }
                    else
                    {
                        int curCharIdx = static_cast<int>(segments[segmentIdx][variantIdx].size()) - 1;
                        int curPatternIdx = static_cast<int>(leaf.second);

                        while (curCharIdx >= 0)
                        {
                            if (curPatternIdx < 0)
                            {
                                const SourceSet newSources = (variantSources & leaf.first);

                                if (newSources.any())
                                    return true;
                            }

                            if (pattern[curPatternIdx] != segments[segmentIdx][variantIdx][curCharIdx])
                                break;

                            curCharIdx -= 1;
                            curPatternIdx -= 1;
                        }

                        if (curPatternIdx < 0)
                        {
                            const SourceSet newSources = (variantSources & leaf.first);

                            if (newSources.any())
                                return true;
                        }

                        if (curCharIdx < 0)
                        {
                            const SourceSet newSources = (variantSources & leaf.first);

                            if (newSources.any())
                            {
                                newLeaves.emplace_back(move(newSources), curPatternIdx);
                            }
                        }
                    }
                }
            }

            leaves = newLeaves;
        }

        segmentIdx -= 1;
    }

    return false;
}

Sopang::SourceSet calcMatchSources(const string *const *segments,
    const int *segmentSizes,
    const Sopang::SourceMap &sourceMap,
    int sourceCount,
    const string &pattern,
    int matchIdx,
    const pair<int, int> &match,
    bool &deterministicSegmentMatch)
{
    using SourceSet = Sopang::SourceSet;
    const int patternCharIdx = static_cast<int>(pattern.size()) - match.second - 2;

    if (patternCharIdx < 0) // The match is fully contained within a single segment.
    {
        if (sourceMap.count(matchIdx) > 0)
        {
            assert(match.first >= 0 and match.first < static_cast<int>(sourceMap.at(matchIdx).size()));
            return sourceMap.at(matchIdx)[match.first];
        }

        deterministicSegmentMatch = true;
        return SourceSet(sourceCount);
    }

    vector<pair<SourceSet, int>> leaves;
    SourceSet rootSources(sourceCount);
    
    if (sourceMap.count(matchIdx) > 0)
    {
        assert(match.first >= 0 and match.first < static_cast<int>(sourceMap.at(matchIdx).size()));
        rootSources = sourceMap.at(matchIdx)[match.first];
    }
    else
    {
        rootSources.set();
    }

    leaves.emplace_back(move(rootSources), patternCharIdx);
    int segmentIdx = static_cast<int>(matchIdx - 1);

    SourceSet res(sourceCount);

    while (not leaves.empty() and segmentIdx >= 0)
    {
        vector<pair<SourceSet, int>> newLeaves;
        newLeaves.reserve(leaves.size() * segmentSizes[segmentIdx]);

        if (segmentSizes[segmentIdx] == 1)
        {
            for (auto &leaf : leaves)
            {
                leaf.second -= segments[segmentIdx][0].size();

                if (leaf.second < 0)
                {
                    res |= leaf.first;
                }
                else
                {
                    newLeaves.push_back(leaf);
                }
            }

            leaves = newLeaves;
        }
        else
        {
            assert(sourceMap.count(segmentIdx) > 0 and sourceMap.at(segmentIdx).size() == static_cast<size_t>(segmentSizes[segmentIdx]));

            for (const auto &leaf : leaves)
            {
                for (int variantIdx = 0; variantIdx < segmentSizes[segmentIdx]; ++variantIdx)
                {
                    const SourceSet &variantSources = sourceMap.at(segmentIdx)[variantIdx];

                    if (segments[segmentIdx][variantIdx].empty())
                    {
                        const SourceSet newSources = (variantSources & leaf.first);

                        if (newSources.any())
                        {
                            newLeaves.emplace_back(move(newSources), leaf.second);
                        }
                    }
                    else
                    {
                        int curCharIdx = static_cast<int>(segments[segmentIdx][variantIdx].size()) - 1;
                        int curPatternIdx = static_cast<int>(leaf.second);

                        while (curCharIdx >= 0)
                        {
                            if (curPatternIdx < 0)
                                break;

                            if (pattern[curPatternIdx] != segments[segmentIdx][variantIdx][curCharIdx])
                                break;

                            curCharIdx -= 1;
                            curPatternIdx -= 1;
                        }

                        if (curPatternIdx < 0)
                        {
                            const SourceSet newSources = (variantSources & leaf.first);

                            if (newSources.any())
                            {
                                res |= newSources;
                                continue;
                            }
                        }

                        if (curCharIdx < 0)
                        {
                            const SourceSet newSources = (variantSources & leaf.first);

                            if (newSources.any())
                            {
                                newLeaves.emplace_back(move(newSources), curPatternIdx);
                            }
                        }
                    }
                }
            }

            leaves = newLeaves;
        }

        segmentIdx -= 1;
    }

    return res;
}

} // namespace (anonymous)

unordered_set<int> Sopang::matchWithSourcesVerify(const string *const *segments,
    int nSegments,
    const int *segmentSizes,
    const Sopang::SourceMap &sourceMap,
    int sourceCount,
    const string &pattern)
{
    const IndexToMatchMap indexToMatch = calcIndexToMatchMap(segments, nSegments, segmentSizes, pattern);
    unordered_set<int> res;

    for (const auto &kv : indexToMatch)
    {
        for (const auto &match : kv.second)
        {
            if (verifyMatch(segments, segmentSizes, sourceMap, sourceCount, pattern, kv.first, match))
            {
                res.insert(kv.first);
                break;
            }
        }
    }

    return res;
}

unordered_map<int, Sopang::SourceSet> Sopang::matchWithSources(const string *const *segments,
    int nSegments,
    const int *segmentSizes,
    const Sopang::SourceMap &sourceMap,
    int sourceCount,
    const string &pattern)
{
    const IndexToMatchMap indexToMatch = calcIndexToMatchMap(segments, nSegments, segmentSizes, pattern);

    unordered_map<int, SourceSet> res;

    for (const auto &kv : indexToMatch)
    {
        for (const auto &match : kv.second)
        {
            bool deterministicSegmentMatch = false;
            const SourceSet curSources = calcMatchSources(segments, segmentSizes, sourceMap, sourceCount, pattern, kv.first, match, deterministicSegmentMatch);

            if (not curSources.empty())
            {
                if (res.count(kv.first) > 0)
                {
                    res.find(kv.first)->second |= curSources;
                }
                else
                {
                    res.emplace(kv.first, move(curSources));
                }
            }
            else if (deterministicSegmentMatch)
            {
                res.emplace(kv.first, sourceCount);
            }
        }
    }

    return res;
}

Sopang::IndexToMatchMap Sopang::calcIndexToMatchMap(const string *const *segments,
    int nSegments,
    const int *segmentSizes,
    const string &pattern)
{
    assert(nSegments > 0 and pattern.size() > 0 and pattern.size() <= wordSize);
    fillPatternMaskBuffer(pattern);

    const uint64_t hitMask = (0x1ULL << (pattern.size() - 1));
    uint64_t D = allOnes;

    // segment index -> [(variant index, char in variant index)]
    IndexToMatchMap res;
    res.reserve(matchMapReserveSize);

    for (int iS = 0; iS < nSegments; ++iS)
    {
        assert(segmentSizes[iS] > 0 and static_cast<size_t>(segmentSizes[iS]) <= dBufferSize);

        for (int iD = 0; iD < segmentSizes[iS]; ++iD)
        {
            dBuffer[iD] = D;

            for (size_t iC = 0; iC < segments[iS][iD].size(); ++iC)
            {
                const char c = segments[iS][iD][iC];

                assert(c > 0 and static_cast<unsigned char>(c) < maskBufferSize);
                assert(alphabet.find(c) != string::npos);

                dBuffer[iD] <<= 1;
                dBuffer[iD] |= maskBuffer[static_cast<unsigned char>(c)];

                if ((dBuffer[iD] & hitMask) == 0x0ULL)
                {
                    res[iS].emplace_back(iD, static_cast<int>(iC));
                }
            }
        }

        D = dBuffer[0];

        for (int iD = 1; iD < segmentSizes[iS]; ++iD)
        {
            D &= dBuffer[iD];
        }
    }

    return res;
}

void Sopang::initCounterPositionMasks()
{
    for (size_t i = 0; i < maxPatternApproxSize; ++i)
    {
        counterPosMasks[i] = (saCounterAllSet << (i * saCounterSize));
    }
}

void Sopang::fillPatternMaskBuffer(const string &pattern)
{
    assert(pattern.size() > 0 && pattern.size() <= wordSize);
    assert(alphabet.size() > 0);

    for (const char c : alphabet)
    {
        assert(c > 0 && static_cast<unsigned char>(c) < maskBufferSize);
        maskBuffer[static_cast<unsigned char>(c)] = allOnes;
    }

    for (size_t iC = 0; iC < pattern.size(); ++iC)
    {
        assert(pattern[iC] > 0 && static_cast<unsigned char>(pattern[iC]) < maskBufferSize);
        maskBuffer[static_cast<unsigned char>(pattern[iC])] &= (~(0x1ULL << iC));
    }
}

void Sopang::fillPatternMaskBufferApprox(const string &pattern)
{
    assert(pattern.size() > 0 and pattern.size() <= wordSize);
    assert(alphabet.size() > 0);

    for (const char c : alphabet)
    {
        assert(c > 0 and static_cast<unsigned char>(c) < maskBufferSize);
        maskBuffer[static_cast<unsigned char>(c)] = 0x0ULL;

        for (size_t iC = 0; iC < pattern.size(); ++iC)
        {
            maskBuffer[static_cast<unsigned char>(c)] |= (0x1ULL << (iC * saCounterSize));
        }
    }

    for (size_t iC = 0; iC < pattern.size(); ++iC)
    {
        assert(pattern[iC] > 0 and static_cast<unsigned char>(pattern[iC]) < maskBufferSize);
        // We zero the bit at the counter position corresponding to the current character in the pattern.
        maskBuffer[static_cast<unsigned char>(pattern[iC])] &= (~(0x1ULL << (iC * saCounterSize)));
    }
}

} // namespace sopang
