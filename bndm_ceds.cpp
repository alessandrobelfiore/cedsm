#include "sopang.hpp"
#include "helpers.hpp"

#include <iostream>
#include <ctime>
#include <unordered_set>

#define SIGMA 256       // Max size of the alphabet

using namespace sopang;

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
