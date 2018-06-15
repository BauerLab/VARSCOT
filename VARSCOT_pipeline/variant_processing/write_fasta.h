// ============================================================================
// Extract variant sequences from an fai index and write them to a fasta file
// ============================================================================
#pragma once

#include <iostream>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

namespace seqan
{
/*!
 * @fn getFastaID
 * @brief Creates a fasta ID for a variant sequence containing information needed to later map sequence back to
 *        genomic regions
 *
 * @signature void getFastaID(fastaID, allVariants, sortedIndex, first, combination, chrName)
 *
 * @param[in,out]   fastaID             Generated fasta ID as string
 * @param[in]       allVariants         Vector of vectors containing all variants from the complete VCF file, variants
 *                                      from the same record are stored in the same vector
 * @param[in]       sortedIndex         Index of position-sorted variants within one chromosome
 * @param[in]       first               First position of range within sortedIndex
 * @param[in]       indexCenterVariant  Variant index marking the variant the range is corresponding to
 * @param[in]       combination         Combination of variants at positions within current range
 * @param[in]       chrName             Real name of the chromosome variants are located
 */
void getFastaID(CharString & fastaID, std::vector<std::vector<VariantSequence> > const & allVariants,
std::vector<unsigned> const & sortedIndex, int first, unsigned const & indexCenterVariant, std::vector<int> combination,
CharString & chrName)
{
    StringSet<CharString> partsOfID;
    reserve(partsOfID, combination.size() * 6 + 5);
    appendValue(partsOfID, chrName);
    appendValue(partsOfID, "_");
    // Position are 0-based results matching final output file of pipeline
    appendValue(partsOfID, (CharString) std::to_string(allVariants[indexCenterVariant][0].start));
    appendValue(partsOfID, "_");

    bool allRef = std::all_of(combination.begin(), combination.end(), [](int i) {return i == -1;});
    if (allRef)
    {
        appendValue(partsOfID, "REF");
    }
    else
    {
        appendValue(partsOfID, "ALT");
        for (unsigned i = 0, j = 1; i < combination.size(); i++, j +=6)
        {
            if (combination[i] != -1)
            {
                appendValue(partsOfID, "_");
                // Position + 1 to get 1-based results matching vcf file
                appendValue(partsOfID, (CharString) std::to_string(allVariants[sortedIndex[first + i]][combination[i]].pos));
                appendValue(partsOfID, "_");
                appendValue(partsOfID, allVariants[sortedIndex[first + i]][combination[i]].ref);
                appendValue(partsOfID, "_");
                appendValue(partsOfID, allVariants[sortedIndex[first + i]][combination[i]].alt);
            }
        }
    }
    fastaID = concat(partsOfID);
}

/*!
 * @fn allCombinations
 * @brief Returns all combinations of variants within a given range (sequences) and fasta IDs
 *
 * @signature void allCombinations(altCombinations, fastaIDs, allVariants, sortedIndex, first, last, chrName)
 *
 * @param[in,out]   altCombinations Vector containing variant (alternative) combinations for one overlapping range
 * @param[in,out]   fastaIDs        Vector of generated fasta IDs as string for combinations of one overlapping range
 * @param[in]       allVariants     Vector of vectors containing all variants from the complete VCF file, variants from
 *                                  the same record are stored in the same vector
 * @param[in]       sortedIndex     Index of position-sorted variants within one chromosome
 * @param[in]       first           First position of range within sortedIndex
 * @param[in]       last            Last position of range within sortedIndex (last not included anymore)
 * @param[in]       indexCenterVariant  Variant index marking the variant the range is corresponding to
 * @param[in]       chrName         Real name of the chromosome where variants are located
 *
 * For each range of variants there might be more than 2 (for 2 alleles) sequence combinations due to unphased variants.
 * As all those variants should be considered every combination of unphased variants (still in the order given by their
 * position and phased variants fixed) needs to be generated. This function returns all combinations of alternatives
 * (as sequence) for a certain range and corresponding fasta IDs.
 */
void allCombinations(std::vector<std::vector<Dna5String> > & altCombinations, std::vector<CharString> & fastaIDs,
std::vector<std::vector<VariantSequence> > const & allVariants, std::vector<unsigned> const & sortedIndex,
unsigned first, unsigned last, unsigned const & indexCenterVariant, CharString & chrName)
{
    unsigned vectorSize = last - first;
    // Stores index in sequence of unphased variants
    std::vector<unsigned> unphasedVariants;
    unphasedVariants.reserve(vectorSize);
    unsigned num = 1;

    // Store the correct variants for both alleles
    std::vector <Dna5String> firstAlleleSequence, secondAlleleSequence;
    firstAlleleSequence.resize(vectorSize);
    secondAlleleSequence.resize(vectorSize);

    // Needed for fasta ID construction
    std::vector<int> indexVariantsFirst, indexVariantsSecond;
    indexVariantsFirst.resize(vectorSize);
    indexVariantsSecond.resize(vectorSize);

    unsigned count = 0;
    // Fill in phased variants and mark unphased variants
    for (unsigned i = first; i < last; i++)
    {
        if (allVariants[sortedIndex[i]][0].allele == -1)
        {
            unphasedVariants.push_back(count);
            num *= 2;
        }
        else
        {
            if (allVariants[sortedIndex[i]].size() == 2)
            {
                firstAlleleSequence[i - first] = allVariants[sortedIndex[i]][0].alt;
                indexVariantsFirst[i - first] = 0;
                secondAlleleSequence[i - first] = allVariants[sortedIndex[i]][1].alt;
                indexVariantsSecond[i - first] = 1;
            }
            else if (allVariants[sortedIndex[i]][0].allele == 0)
            {
                firstAlleleSequence[i - first] = allVariants[sortedIndex[i]][0].alt;
                indexVariantsFirst[i - first] = 0;
                secondAlleleSequence[i - first] = allVariants[sortedIndex[i]][0].ref;
                indexVariantsSecond[i - first] = -1;
            }
            else if (allVariants[sortedIndex[i]][0].allele == 1)
            {
                firstAlleleSequence[i - first] = allVariants[sortedIndex[i]][0].ref;
                indexVariantsFirst[i - first] = -1;
                secondAlleleSequence[i - first] = allVariants[sortedIndex[i]][0].alt;
                indexVariantsSecond[i - first] = 0;
            }
            else
            {
                firstAlleleSequence[i - first] = allVariants[sortedIndex[i]][0].alt;
                indexVariantsFirst[i - first] = 0;
                secondAlleleSequence[i - first] = allVariants[sortedIndex[i]][0].alt;
                indexVariantsFirst[i - first] = 0;
            }
        }
        count++;
    }

    altCombinations.reserve(2 * num);
    fastaIDs.reserve(2 * num);
    CharString fastaID;

    if (unphasedVariants.size() > 0)
    {
        // Stack keeps track of indices in all vectors from the range (unphased variants)
        std::vector<int> stack = {-1};
        while (!stack.empty())
        {
            ++stack.back(); // Last index + 1
            if (stack.back() >= 2)
            {
                // Last index in stack is too large, index before must be increased by one
                stack.pop_back();
            }
            else if (stack.size() < unphasedVariants.size())
            {
                // Stack is not containing a valid tuple yet, so another subindex need to be added
                stack.push_back(-1);
            }
            else
            {
                // Save new tuple
                for (unsigned i = 0; i < stack.size(); ++i)
                {
                    if (allVariants[sortedIndex[first + unphasedVariants[i]]].size() == 2)
                    {
                        firstAlleleSequence[unphasedVariants[i]] = allVariants[sortedIndex[first + unphasedVariants[i]]][stack[i]].alt;
                        indexVariantsFirst[unphasedVariants[i]] = stack[i];
                        secondAlleleSequence[unphasedVariants[i]] = allVariants[sortedIndex[first + unphasedVariants[i]]][stack[i]].alt;
                        indexVariantsSecond[unphasedVariants[i]] = stack[i];
                    }
                    else // Assuming 0/0, 1/1, 2/2 are not possible because excluded in processVCF from "unphased" (we don't need the combinations)
                    {
                        if (stack[i] == 0)
                        {
                            firstAlleleSequence[unphasedVariants[i]] = allVariants[sortedIndex[first + unphasedVariants[i]]][0].ref;
                            indexVariantsFirst[unphasedVariants[i]] = -1; // Indicates ref for fasta IDs
                            secondAlleleSequence[unphasedVariants[i]] = allVariants[sortedIndex[first + unphasedVariants[i]]][0].ref;
                            indexVariantsSecond[unphasedVariants[i]] = -1;
                        }
                        else
                        {
                            firstAlleleSequence[unphasedVariants[i]] = allVariants[sortedIndex[first + unphasedVariants[i]]][0].alt;
                            indexVariantsFirst[unphasedVariants[i]] = 0;
                            secondAlleleSequence[unphasedVariants[i]] = allVariants[sortedIndex[first + unphasedVariants[i]]][0].alt;
                            indexVariantsSecond[unphasedVariants[i]] = 0;
                        }
                    }
                }
                altCombinations.push_back(firstAlleleSequence);
                getFastaID(fastaID, allVariants, sortedIndex, first, indexCenterVariant, indexVariantsFirst, chrName);
                fastaIDs.push_back(fastaID);

                if (indexVariantsFirst != indexVariantsSecond)
                {
                    altCombinations.push_back(secondAlleleSequence);
                    getFastaID(fastaID, allVariants, sortedIndex, first, indexCenterVariant, indexVariantsSecond, chrName);
                    fastaIDs.push_back(fastaID);
                }
            }
        }
    }
    else
    {
        altCombinations.push_back(firstAlleleSequence);
        getFastaID(fastaID, allVariants, sortedIndex, first, indexCenterVariant, indexVariantsFirst, chrName);
        fastaIDs.push_back(fastaID);

        if (indexVariantsFirst != indexVariantsSecond)
        {
            altCombinations.push_back(secondAlleleSequence);

            getFastaID(fastaID, allVariants, sortedIndex, first, indexCenterVariant, indexVariantsSecond, chrName);
            fastaIDs.push_back(fastaID);
        }
    }
}

/*!
 * @fn extractSequenceFromIndex
 * @brief Extracts a sequence from an fai index based on chromosome, start and end.
 *
 * @signature void extractSequenceFromIndex(sequence, faiIndex, id, start, end)
 *
 * @param[in,out]   sequence        Extracted sequence
 * @param[in]       faiIndex        The fai index to be searched
 * @param[in]       id              The fasta ID were the sequence is located
 * @param[in]       start           Start position of the sequence (0-based)
 * @param[in]       end             End position of the sequence (not included, half open intervals)
 *
 * @throw Exception if sequence cannot be read in
 */
void extractSequenceFromIndex(Dna5String & sequence, FaiIndex & faiIndex, CharString & id, unsigned start, unsigned end)
{
    // Translate sequence name to index.
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, id))
    {
        throw std::out_of_range("ERROR: Index out of range.");
    }

    // Make sure start <= end <= sequenceLength
    if (start > sequenceLength(faiIndex, idx))
        start = sequenceLength(faiIndex, idx);
    if (end > sequenceLength(faiIndex, idx))
        end = sequenceLength(faiIndex, idx);
    if (start > end)
        end = start;

    try
    {
        readRegion(sequence, faiIndex, idx, start, end);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }
}

/*!
 * @fn generateVariantSequences
 * @brief Generates all fasta sequences of variants for one range of one chromosome
 *
 * @signature void generateVariantSequences(fastaSequences, fastaIDs, faiIndex, allVariants, sortedIndex, chrName,
 *                                          range, seqLength)
 *
 * @param[in,out]   fastaSequences      All fasta sequences of variants for one range of one chromosome
 * @param[in,out]   fastaIDs            Corresponding fasta IDs
 * @param[in]       faiIndex            The fai index to be searched
 * @param[in]       allVariants         Vector of vectors containing all variants from the complete VCF file, variants
 *                                      from the same record are stored in the same vector
 * @param[in]       sortedIndex         Index of position-sorted variants within one chromosome
 * @param[in]       chrName             Real name of the chromosome where variants are located
 * @param[in]       range               Current range processed
 * @param[in]       indexCenterVariant  Variant index marking the variant the range is corresponding to
 * @param[in]       seqLength           Length of the target sequence that should be searched later (normally 23 bp)
 *
 * This function generates sequences which contain one combination of variants of one range. The sequence has a start
 * position such that a window of length seqLength can be slidden over the sequence in a way that the left-most variant
 * starts as last base of the first window. This works analogous for the end position. This might cause duplicates in
 * both directions:
 * - In case of a deletion at the first variant the first seqLength bp are actually reference sequence because the left-
 * most base in an indel VCF record is always reference.
 * - In case of an indel that has reference base at first and last position of alternative the same happens to the end
 * of the sequence
 * As those cases might be rare this is accepted and later in the pipeline duplicates are removed.
 * The function first extracts reference sequence needed in the begin and end (if not variant is at the ends) as well
 * as in between of variants and afterwards inserts aternative sequences and creates the whole sequence.
*/
void generateVariantSequences(std::vector<Dna5String> & fastaSequences, std::vector<CharString> & fastaIDs,
FaiIndex & faiIndex, std::vector<std::vector<VariantSequence> > const & allVariants,
std::vector<unsigned> & sortedIndex, CharString & chrName, Pair<unsigned, unsigned> & range,
unsigned const & indexCenterVariant, unsigned const & seqLength)
{
    // Size of the stringset depends on whether variants lie on the edges or not
    // First check if first and/or last variant lie on edges and need to be considered

    bool startVariant = false;
    bool endVariant = false;

    // Start and first variant
    if (allVariants[indexCenterVariant][0].start > allVariants[sortedIndex[getValueI1(range)]][0].pos)
    {
        startVariant = true;
    }

    if (allVariants[indexCenterVariant][0].end == allVariants[sortedIndex[getValueI2(range) - 1]][0].pos)
    {
        endVariant = true;
    }

    // Define StringSet for sequences
    StringSet<Dna5String> sequenceBase;
    // Set parameters for for loop
    // If variant is at first or last position than no reference is needed to the left/right
    unsigned refSeqStart, rangeStart, rangeEnd;


    if (startVariant && endVariant)
    {
        resize(sequenceBase, 2 * (getValueI2(range) - getValueI1(range)) - 1);
        refSeqStart = 1;
        rangeStart = getValueI1(range) + 1;
        rangeEnd = getValueI2(range); // Half open interval
    }
    else if (startVariant)
    {
        resize(sequenceBase, 2 * (getValueI2(range) - getValueI1(range)));
        refSeqStart = 1;
        rangeStart = getValueI1(range) + 1;
        rangeEnd = getValueI2(range) + 1; // Half open interval
    }
    else if (endVariant)
    {
        resize(sequenceBase, 2 * (getValueI2(range) - getValueI1(range)));
        refSeqStart = 0;
        rangeStart = getValueI1(range);
        rangeEnd = getValueI2(range); // Half open interval
    }
    else
    {
        resize(sequenceBase, 2 * (getValueI2(range) - getValueI1(range)) + 1);
        refSeqStart = 0;
        rangeStart = getValueI1(range);
        rangeEnd = getValueI2(range) + 1; // Half open interval
    }

    // First fill in the reference sequences
    unsigned beginSeq, endSeq;

    // Different sequence lengths depending on first, last or in-between variant
    // Need to take into account that sequence might start or end in the middle of a variant
    for (unsigned i = rangeStart, j = refSeqStart; i < rangeEnd; i++, j += 2)
    {
        if (j == 0)
        {
            beginSeq = allVariants[indexCenterVariant][0].start;
            endSeq = allVariants[sortedIndex[i]][0].pos;
        }
        else if (i == getValueI2(range))
        {
            beginSeq = allVariants[sortedIndex[i - 1]][0].pos + length(allVariants[sortedIndex[i - 1]][0].ref);
            endSeq = allVariants[indexCenterVariant][0].end;
        }
        else
        {
            beginSeq = allVariants[sortedIndex[i - 1]][0].pos + length(allVariants[sortedIndex[i - 1]][0].ref);
            endSeq = allVariants[sortedIndex[i]][0].pos;
        }
        extractSequenceFromIndex(sequenceBase[j], faiIndex, chrName, beginSeq, endSeq);
    }

    // Then generate every combination of variant sequences
    std::vector<std::vector<Dna5String> > altCombinations;
    allCombinations(altCombinations, fastaIDs, allVariants, sortedIndex, getValueI1(range), getValueI2(range),
                    indexCenterVariant, chrName);
    fastaSequences.resize(altCombinations.size());
    for (unsigned i = 0; i < altCombinations.size(); i++)
    {
        for (unsigned j = 0, k = (1u - refSeqStart); j < altCombinations[i].size(); j++, k += 2)
        {
            assignValue(sequenceBase, k, altCombinations[i][j]);
        }
        fastaSequences[i] = concat(sequenceBase);
    }
}

/*!
 * @fn writeFastaFile
 * @brief Writes variant sequences which are created on the fly to a fasta file
 *
 * @signature void writeFastaFile(outputPath, indexPath, allVariants, sortedIndexAllChr, allOverlapRegions, chrTable,
 *                                seqLength)
 *
 * @param[in]       outputPath              Output path for fasta file
 * @param[in]       indexPath               Path for the genome index should be built of or index exists in same
 *                                          directory
 * @param[in]       allVariants             Vector of vectors containing all variants from the complete VCF file,
 *                                          variants
 *                                          from the same record are stored in the same vector
 * @param[in]       sortedIndexAllChr       Index of position-sorted variants for all chromosomes (reading order)
 * @param[in]       allOverlapRegions       Vector containing half open intervals marking the largest overlap of
 *                                          subsequent variants for all chromosomes (reading order)
 * @param[in]       allIndexCenterVariants  Vector containing for each range the variant index from sorted index from
 *                                          which range was computed for all chromosomes (reading order)
 * @param[in]       chrTable                Look up table containing the chromosome names at the corresponding index
 * @param[in]       seqLength               Length of the target sequence that should be searched later
 *                                          (normally 23 bp)
*/
void writeFastaFile(char const * outputPath, char const * indexPath,
std::vector<std::vector<VariantSequence> > const & allVariants, std::vector<std::vector<unsigned> > & sortedIndexAllChr,
std::vector<std::vector<Pair<unsigned, unsigned> > > & allOverlapRegions,
std::vector<std::vector<unsigned> > const & allIndexCenterVariants, std::vector<CharString> & chrTable,
unsigned const & seqLength)
{
    SeqFileOut seqOut;

    if (!open(seqOut, outputPath))
    {
        throw std::runtime_error("ERROR: Could not open single FASTA output file.");
    }

    // Try to load index and create on the fly if necessary
    FaiIndex faiIndex;
    if (!open(faiIndex, indexPath))
    {
        if (!build(faiIndex, indexPath))
        {
            throw std::runtime_error("ERROR: Index could not be loaded or built.");
        }
        if (!save(faiIndex)) // Name is stored from when reading.
        {
            throw std::runtime_error("ERROR: Index could not be written do disk.");
        }
    }

    // For each chromosome, process ranges and create overlapping sequences on the fly
    try
    {
        for (unsigned i = 0; i < allOverlapRegions.size(); i++)
        {
            for (unsigned j = 0; j < allOverlapRegions[i].size(); j++)
            {
                std::vector<Dna5String> fastaSequences;
                std::vector<CharString> fastaIDs;
                generateVariantSequences(fastaSequences, fastaIDs, faiIndex, allVariants, sortedIndexAllChr[i],
                                         chrTable[i], allOverlapRegions[i][j], allIndexCenterVariants[i][j], seqLength);
                writeRecords(seqOut, fastaIDs, fastaSequences);
            }
        }
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return;
    }
}

}
