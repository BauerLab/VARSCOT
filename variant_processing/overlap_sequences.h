// ============================================================================
// Compute overlap sequences of subsequent variants
// ============================================================================
#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include <seqan/parallel.h>

namespace seqan
{
/*!
 * @fn findMaxOverlap
 * @brief Find maximum ranges of overlapping subsequent variants within one chromosome
 *
 * @signature findMaxOverlap(overlapRegions, allVariants, sortedIndex, seqLength)
 *
 * @param[in,out]   overlapRegions      Vector containing half open intervals marking the largest overlap of subsequent
 *                                      variants within one chromosome
 * @param[in,out]   indexCenterVariants Vector containing for each range the variant index from sorted index from which
 *                                      range was computed
 * @param[in]       allVariants         Vector of vectors containing all variants from the complete VCF file, variants
 *                                      from the same record are stored in the same vector
 * @param[in]       sortedIndex         Index of position-sorted variants within one chromosome
 * @param[in]       seqLength           Length of the target sequence that should be searched later (normally 23 bp)
 *
 * This function computes ranges variants that overlap and lie within a region of 22 bp +/- every variant.
 * Region is extended for largest deletion at all considered positions to make sure no variant is missing.
 * This produces duplicates as windows are probably overlapping and if variants exist where both left-most and right-most
 * bases in the alternative sequence belong to the reference. As this is a rare case it is accepted and later
 * in the pipeline duplicates are filtered out.
 */
void findMaxOverlap(std::vector<Pair<unsigned, unsigned> > & overlapRegions, std::vector<unsigned> & indexCenterVariants,
std::vector<std::vector<VariantSequence> > & allVariants, std::vector<unsigned> const & sortedIndex,
unsigned const & seqLength)
{
    // Check for largest deletion at each position
    // Window is increased by largest deletion at all considered variant positions
    std::vector<unsigned> maxDeletion(sortedIndex.size(), 0);
    for (unsigned i = 0; i < sortedIndex.size(); i++)
    {
        for (auto j : allVariants[sortedIndex[i]])
        {
            if (j.variantType == 2)
            {
                if ((length(j.ref) - length(j.alt)) > maxDeletion[i])
                    maxDeletion[i] = length(j.ref) - length(j.alt);
            }
        }
    }

    overlapRegions.reserve(sortedIndex.size());
    indexCenterVariants.reserve(sortedIndex.size());
    Pair<unsigned, unsigned> range;

    setValueI1(range, 0);
    setValueI2(range, 0);
    unsigned windowSizeLeft, windowSizeRight;

    for (unsigned i = 0; i < sortedIndex.size(); i++)
    {
        int indexLeft, indexRight;

        // First check whether left check is necessary or not
        // If the preceeding variant did not include the current than don't need to check left
        if (getValueI2(range) > i)
        {
            // First check right to avoid duplicates
            // If no addtional variant lies within right region than this variants' whole region is already included in preceeding range
            indexRight = getValueI2(range);

            // Need to include all earlier maximal deletions into window
            // TODO check if this can be saved from the earlier variants
            // Only need to consider maximum length of insertion or deletion
            windowSizeRight = seqLength + maxDeletion[i];
            if (indexRight < sortedIndex.size())
            {
                for (unsigned del = i + 1; del <= indexRight; del++)
                {
                    windowSizeRight += maxDeletion[del];
                }
            }

            while (indexRight < sortedIndex.size() &&
            (allVariants[sortedIndex[indexRight]][0].pos - allVariants[sortedIndex[i]][0].pos) < windowSizeRight)
            {
                if (indexRight < sortedIndex.size()) // TODO Try to get rid of this check
                {
                    windowSizeRight += maxDeletion[indexRight];
                }
                indexRight++;
            }
            // If no new variant is found edit end position of preceeding range center variant
            // Can be done without issues as no variant was found within window size and therefore there must be only reference
            if (indexRight == getValueI2(range))
            {
                for (unsigned j = 0; j < allVariants[indexCenterVariants.back()]. size(); j++)
                {
                    allVariants[indexCenterVariants.back()][j].end = allVariants[sortedIndex[i]][j].pos + windowSizeRight;
                }
                continue;
            }
            setValueI2(range, indexRight); // half open intervals, indexRight is position of the first variant that lies not within range anymore

            // Left check needed
            indexLeft = i - 1;
            windowSizeLeft = seqLength + maxDeletion[indexLeft];

            while (indexLeft >= 0 &&
                (allVariants[sortedIndex[i]][0].pos - allVariants[sortedIndex[indexLeft]][0].pos) < windowSizeLeft)
            {
                indexLeft--;
                windowSizeLeft += maxDeletion[indexLeft];
            }

            // Check if preceeding range is included in current range
            // If so, than change preceeding range (new end position) and continue (don't save a new range)
            if ((indexLeft + 1) == getValueI1(range))
            {
                for (unsigned j = 0; j < allVariants[indexCenterVariants.back()]. size(); j++)
                {
                    allVariants[indexCenterVariants.back()][j].end = allVariants[sortedIndex[i]][j].pos + windowSizeRight;
                }
                setValueI2(overlapRegions[overlapRegions.size() - 1], indexRight);
                continue;
            }
            setValueI1(range, indexLeft + 1); // Because while loop terminates such that indexLeft points at the variant that is already out of range
        }
        else
        {
            // Check right
            windowSizeRight = seqLength + maxDeletion[i];
            indexRight = i + 1;

            while (indexRight < sortedIndex.size() &&
                (allVariants[sortedIndex[indexRight]][0].pos - allVariants[sortedIndex[i]][0].pos) < windowSizeRight)
            {
                if (indexRight < sortedIndex.size()) // TODO Try to get rid of this check
                {
                    windowSizeRight += maxDeletion[indexRight];
                }
                indexRight++;
            }
            setValueI2(range, indexRight); // half open intervals, indexRight is position of the first variant that lies not within range anymore

            // Only reference and no variants to the left
            windowSizeLeft = seqLength;
            setValueI1(range, i);
        }
        overlapRegions.push_back(range);
        indexCenterVariants.push_back(sortedIndex[i]);

        // Set end position of sequences
        for (unsigned j = 0; j < allVariants[sortedIndex[i]]. size(); j++)
        {
            allVariants[sortedIndex[i]][j].start = allVariants[sortedIndex[i]][j].pos - windowSizeLeft + 1;
            allVariants[sortedIndex[i]][j].end = allVariants[sortedIndex[i]][j].pos + windowSizeRight;
        }
    }
}

/*!
 * @fn getVariantOverlapRanges
 * @brief Find maximum ranges of overlapping subsequent variants for all chromosomes and variants
 *
 * @signature void getVariantOverlapRanges(allOverlapRegions, sortedIndexAllChr, chrNumber, allVariants, seqLength,
 *                                         threads)
 *
 * @param[in,out]   allOverlapRegions       Vector containing half open intervals marking the largest overlap of
 *                                          subsequent variants for all chromosomes (reading order)
 * @param[in,out]   allIndexCenterVariants  Vector containing for each range the variant index from sorted index from
 *                                          which range was computed for all chromosomes (reading order)
 * @param[in,out]   sortedIndexAllChr       Index of position-sorted variants for all chromosomes (reading order)
 * @param[in]       chrNumber               Number of chromosomes read from VCF file (might not be all chromosomes)
 * @param[in]       allVariants             Vector of vectors containing all variants from the complete VCF file,
 *                                          variants from the same record are stored in the same vector
 * @param[in]       seqLength               Length of the target sequence that should be searched later
 *                                          (normally 23 bp)
 * @param[in]       threads                 Number of threads used for parallelization
 */
void getVariantOverlapRanges(std::vector<std::vector<Pair<unsigned, unsigned> > > & allOverlapRegions,
std::vector<std::vector<unsigned> > & allIndexCenterVariants, std::vector<std::vector<unsigned> > & sortedIndexAllChr,
unsigned const & chrNumber, std::vector<std::vector<VariantSequence> > & allVariants, unsigned const & seqLength, unsigned threads)
{
    // Get indices by chromosome:
    // First, count number of occurrences of each chromosome - O(n) if n is length of input vector
    // Reserve space                                         - O(m) if m is number of chromosomes,
    //                                                         at most 24 (chr 1-22, X and Y)
    // Second, get indices and insert with push_back         - O(n) and never needs to request new space
    //                                                         because already reserved

    // Two vectors:
    // One containing indices by chromosome (index is equal to chromosome index from cache)
    // One containing number of samples per chromosome for space reservation

    sortedIndexAllChr.resize(chrNumber);

    std::vector<unsigned> countIndices(chrNumber, 0);

    // Count number of samples per chromosome
    for (unsigned i = 0; i < allVariants.size(); i++)
    {
        countIndices[allVariants[i][0].chr]++;
    }

    // Reserve space
    for (unsigned i = 0; i < countIndices.size(); i++)
    {
        sortedIndexAllChr[i].reserve(countIndices[i]);
    }

    // Insert indices
    for (unsigned i = 0; i < allVariants.size(); i++)
    {
        sortedIndexAllChr[allVariants[i][0].chr].push_back(i);
    }

    // Sort indices by start
    // Find all windows where subsequent variants lie within a 23 bp region
    // Process all chromosomes in parallel

    allOverlapRegions.resize(chrNumber);
    allIndexCenterVariants.resize(chrNumber);

    omp_set_num_threads(threads);

    #pragma omp parallel for schedule(static)
    for (unsigned i = 0; i < chrNumber; i++)
    {
        printf("thread id: %i\n", omp_get_thread_num());
        std::sort(sortedIndexAllChr[i].begin(), sortedIndexAllChr[i].end(),
            [&](unsigned v1, unsigned v2)
            {
                return allVariants[v1][0].pos < allVariants[v2][0].pos;
            });
        findMaxOverlap(allOverlapRegions[i], allIndexCenterVariants[i], allVariants, sortedIndexAllChr[i], seqLength);
    }
}

}
