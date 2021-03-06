#ifndef SOPANG_HPP
#define SOPANG_HPP

#include "bitset.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef SOPANG_WHITEBOX
#define SOPANG_WHITEBOX
#endif

namespace sopang
{

class Sopang
{
private:
    /** Maximum number of sources (upper bound on source set size). */
    static constexpr int maxSourceCount = 5'120;

public:
    using SourceSet = BitSet<maxSourceCount>;
    using SourceMap = std::unordered_map<int, std::vector<SourceSet>>;

    Sopang(const std::string &alphabet, int sourceCount);
    ~Sopang();

    void shiftSources(Sopang::SourceSet* sources [], int size);

    std::unordered_set<int> match(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const std::string &pattern);

    std::unordered_set<int> matchCEDS(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const SourceMap &sourceMap,
        int sourceCount,
        const std::string &pattern);

    std::unordered_set<int> matchBNDM_CEDS(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const SourceMap &sourceMap,
        int sourceCount,
        const std::string &pattern);

/*     std::unordered_set<int> matchCEDS2(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const SourceMap &sourceMap,
        int sourceCount,
        const std::string &pattern); */

    std::unordered_set<int> matchApprox(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const std::string &pattern,
        int k);

    std::unordered_set<int> matchWithSourcesVerify(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const SourceMap &sourceMap,
        int sourceCount,
        const std::string &pattern);

    std::unordered_map<int, SourceSet> matchWithSources(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const SourceMap &sourceMap,
        int sourceCount,
        const std::string &pattern);

private:
    using IndexToMatchMap = std::unordered_map<int, std::vector<std::pair<int, int>>>;

    IndexToMatchMap calcIndexToMatchMap(const std::string *const *segments,
        int nSegments,
        const int *segmentSizes,
        const std::string &pattern);

    void initCounterPositionMasks();

    void initIndexHitMask();
    void initSourcesVector(int sourceCount);

    void fillPatternMaskBuffer(const std::string &pattern);
    void fillPatternMaskBufferApprox(const std::string &pattern);

    /** Buffer size for processing segment variants, the size of the largest segment (i.e. the number of variants)
     * from the input file cannot be larger than this value. */
    static constexpr size_t dBufferSize = 5'120;
    /** Buffer size for Shift-Or masks for the input alphabet, must be larger than the largest input character ASCII code, 
     * up to 'Z' = 90. */
    static constexpr size_t maskBufferSize = 91;
    /** Word size (in bits) used by the Shift-Or algorithm. */
    static constexpr size_t wordSize = 64;

    /** Maximum pattern size for approximate search. */
    static constexpr size_t maxPatternApproxSize = 12;
    /** Shift-Add counter size in bits. */
    static constexpr size_t saCounterSize = 5;
    
    /** Full single Shift-Add counter indicating no match. */
    static constexpr uint64_t saFullCounter = 0x10ULL;
    /** Single Shift-Add counter with all bits set. */
    static constexpr uint64_t saCounterAllSet = 0x20ULL - 0x1ULL;
    /** Counter with all bits set. */
    static constexpr uint64_t allOnes = ~(0x0ULL);

    /** Initial memory reserve size for a map storing matches for verification with sources. */
    static constexpr size_t matchMapReserveSize = 32;

    uint64_t counterPosMasks[maxPatternApproxSize];

    uint64_t *dBuffer;
    uint64_t maskBuffer[maskBufferSize];

    /* CEDS */
    uint64_t indexHitMask[wordSize];
    Sopang::SourceSet* S[wordSize];
    Sopang::SourceSet* sBuffer[maxSourceCount][wordSize];

    const std::string alphabet;
    int sourceCount;

    SOPANG_WHITEBOX
};

} // namespace sopang

#endif // SOPANG_HPP
