// ============================================================================
// Dunctions to filter bam output from mapping against reference and SNP genome
// ============================================================================

#pragma once

#include <string>
#include <vector>

#include <seqan/parallel.h>
#include <seqan/seq_io.h>

#include "extract_fasta_ontargets.h"

namespace seqan
{
/*!
 * @struct PotentialOffTarget
 * @brief Stores information needed to define a potential off-target
 *
 * @signature struct PotentialOffTarget;
 *
 * The PotentialOffTarget stores information about a potential off-target including the corresponding on-target,
 * the chromosome and mapping position, strand, sequence, positions of mismatches (if existing) and the type of snps
 * included in that sequence (REF, SUB, INS, DEL).
 * Note: Positions are 0-based (internal).
 */
struct PotentialOffTarget
{
    CharString chr;
    CharString target;
    CharString snpType;
    Dna5String sequence;
    std::vector<int> mismatchPos;
    unsigned pos;
    char strand;
};

//!@brief Comparison operator for struct
inline bool comp(PotentialOffTarget const & pot1, PotentialOffTarget const & pot2)
{
    return pot1.target == pot2.target &&
           pot1.chr == pot2.chr &&
           pot1.pos == pot2.pos &&
           pot1.strand == pot2.strand &&
           pot1.sequence == pot2.sequence &&
           pot1.mismatchPos == pot2.mismatchPos &&
           pot1.snpType == pot2.snpType;
}

/*!
 * @fn filterRefAlignment
 * @brief Filter out all mappings to the reference that map to SNP regions (by that also filters out duplicates) or have
 * not a NGG PAM
 *
 * @signature void filterRefAlignment(validRefAlignment, sortedIndexAllChr, chrMap, snpInfoTable, offTargets, seqLength,
 *                                    threads)
 *
 * @param[in,out]   validRefAlignment   Index of valid potential off-targets not mapping into SNP regions
 * @param[in]       sortedIndexAllChr   Sorted index of snpInfoTable by chromosome and position
 * @param[in]       chrMap              Set containing all chromosomes of SNP regions and the corresponding index in
 *                                      sortedIndexAllChr
 * @param[in]       snpInfoTable        Contains all parts of the fasta ID of SNP sequences and length of corresponding
 *                                      fasta sequence
 * @param[in]       offTargets          Potential off-target list
 * @param[in]       onTargets           On-target list
 * @param[in]       seqLength           Length of the target sequences that were searched
 * @param[in]       threads             Number of threads used for parallelization
 */
void filterRefAlignment(std::vector<unsigned> & validRefAlignment,
std::vector<std::vector<unsigned> > const & sortedIndexAllChr, std::map<CharString, unsigned> const & chrMap,
std::vector<StringSet<CharString> > const & snpInfoTable, std::vector<PotentialOffTarget> const & offTargets,
std::map<CharString, PotentialOffTarget> const & onTargets,unsigned const & seqLength, unsigned threads)
{
    validRefAlignment.reserve(offTargets.size());

    // Filter out all reference mappings that lie within SNP regions or are on-targets
    std::vector<unsigned> indexValid(offTargets.size(), 0);
    omp_set_num_threads(threads);

    #pragma omp parallel for schedule(static)
    for (unsigned i = 0; i < offTargets.size(); i++)
    {
        bool valid = true;
        // Check if on-target
        if (comp(offTargets[i], onTargets.at(offTargets[i].target)))
        {
            valid = false;
        }

        // Check location
        // To be removed chromosome must be the same and start and end must lie within a SNP region
        // Every mapping that only overlaps in parts with a SNP regions is valid
        if (valid)
        {
            std::map<CharString, unsigned>::const_iterator it = chrMap.find(offTargets[i].chr);
            if (it != chrMap.end())
            {
                for (unsigned j = 0; j < sortedIndexAllChr[chrMap.at(offTargets[i].chr)].size(); j++)
                {
                    StringSet<CharString> const & currentInfo = snpInfoTable[sortedIndexAllChr[chrMap.at(offTargets[i].chr)][j]];
                    if (offTargets[i].pos >= atoi(toCString(currentInfo[1])) &&
                       (offTargets[i].pos + seqLength) <= (atoi(toCString(currentInfo[1]))) + atoi(toCString(currentInfo[length(currentInfo) - 1])))
                    {
                        valid = false;
                        break;
                    }
                }
            }
        }
        if (valid)
        {
            indexValid[i] = 1;
        }
    }

    for (unsigned i = 0; i < indexValid.size(); i++)
    {
        if (indexValid[i] == 1)
        {
            validRefAlignment.push_back(i);
        }
    }
}

/*!
 * @fn sortSnpRegionsByChr
 * @brief Sort SNP regions by chromosome and position
 *
 * @signature void sortSnpRegionsByChr(sortedIndexAllChr, chrMap, snpInfoTable, threads)
 *
 * @param[in,out]   sortedIndexAllChr   Sorted index of snpInfoTable by chromosome and position
 * @param[in,out]   chrMap              Set containing all chromosomes of SNP regions and the corresponding index in
 *                                      sortedIndexAllChr
 * @param[in]       snpInfoTable        Contains all parts of the fasta ID of SNP sequences and length of corresponding
 *                                      fasta sequence
 * @param[in]       threads             Number of threads used for parallelization
 */
void sortSnpRegionsByChr(std::vector<std::vector<unsigned> > & sortedIndexAllChr,
std::map<CharString, unsigned> const & chrMap, std::vector<StringSet<CharString> > const & snpInfoTable, unsigned threads)
{
    sortedIndexAllChr.resize(chrMap.size()); // All autosomes and sex chromosomes

    std::vector<unsigned> countIndices(chrMap.size(), 0);

    // Count number of samples per chromosome
    for (unsigned i = 0; i < snpInfoTable.size(); i++)
    {
        countIndices[chrMap.at(snpInfoTable[i][0])]++;
    }

    // Reserve space
    for (unsigned i = 0; i < countIndices.size(); i++)
    {
        sortedIndexAllChr[i].reserve(countIndices[i]);
    }

    // Insert indices
    for (unsigned i = 0; i < snpInfoTable.size(); i++)
    {
        sortedIndexAllChr[chrMap.at(snpInfoTable[i][0])].push_back(i);
    }

    // Sort indices by start
    omp_set_num_threads(threads);

    #pragma omp parallel for schedule(static)
    for (unsigned i = 0; i < chrMap.size(); i++)
    {
        std::sort(sortedIndexAllChr[i].begin(), sortedIndexAllChr[i].end(),
            [&](unsigned v1, unsigned v2)
            {
                return atoi(toCString(snpInfoTable[v1][1])) < atoi(toCString(snpInfoTable[v2][1]));
            });
    }
}

/*!
 * @fn getSnpType
 * @brief Get type of SNPs in sequence from fasta ID
 *
 * @signature void getSnpType(snpType, fastaID)
 *
 * @param[in,out]   snpType             Type of SNPs included in a fasta sequence (from fasta ID)
 * @param[in]       fastaID             FastaID parts in StringSet
 * @param[in]       offTargetPos        Starting position of potential off-target
 * @param[in]       seqLength           Length of the target sequences that were searched
 */
void getSnpType(CharString & snpType, StringSet<CharString> const & fastaID, unsigned const & offTargetPos,
unsigned const & seqLength)
{
    StringSet<CharString> variants;
    reserve(variants, (length(fastaID) - 3) / 3);

    appendValue(variants, "VAR_");
    appendValue(variants, fastaID[0]);
    appendValue(variants, "_");
    for (unsigned i = 3; i < length(fastaID); i += 3)
    {
        // Substitution
        if (length(fastaID[i + 1]) == length(fastaID[i + 2]))
        {
            // Check if start position of variant lies within off-target
            // Positions in fasta header are 0-based
            if (offTargetPos <= atoi(toCString(fastaID[i])) &&
               (offTargetPos + seqLength) > atoi(toCString(fastaID[i])))
            {
                appendValue(variants, fastaID[i]);
                appendValue(variants, ',');
            }
        }
        // Insertion
        else if (length(fastaID[i + 1]) < length(fastaID[i + 2]))
        {
            // Check if either first position + 1 (first base is reference) or last position of alternative
            // lies within off-target
            if ((offTargetPos <= (atoi(toCString(fastaID[i])) + 1) &&
                (offTargetPos + seqLength) > (atoi(toCString(fastaID[i])) + 1)) ||
                (offTargetPos <= (atoi(toCString(fastaID[i])) + length(fastaID[i + 2]) - 1) &&
                (offTargetPos + seqLength) > (atoi(toCString(fastaID[i])) + length(fastaID[i + 2]) - 1)))
            {
                appendValue(variants, fastaID[i]);
                appendValue(variants, ',');
            }
        }
        // Deletion
        else
        {
            // Check if either first position + 1 (first base is reference) or last position of reference
            // lies within off-target
            if ((offTargetPos <= (atoi(toCString(fastaID[i])) + 1) &&
                (offTargetPos + seqLength) > (atoi(toCString(fastaID[i])) + 1)) ||
                (offTargetPos <= (atoi(toCString(fastaID[i])) + length(fastaID[i + 1]) - 1) &&
                (offTargetPos + seqLength) > (atoi(toCString(fastaID[i])) + length(fastaID[i + 1]) - 1)))
            {
                appendValue(variants, fastaID[i]);
                appendValue(variants, ',');
            }
        }
    }

    if (length(variants) > 3)
    {
        CharString combinedVariants = concat(variants);
        snpType = prefix(concat(variants), length(combinedVariants) - 1);
    }
}

/*!
 * @fn filterSnpAlignment
 * @brief Filter out all mappings to the SNP genome that map to N regions
 *
 * @signature void filterRefAlignment(validSnpAlignment, offTargets, onTargets, snpInfoTable, seqLength)
 *
 * @param[in,out]   validRefAlignment   Index of valid potential off-targets not mapping into SNP regions
 * @param[in,out]   offTargets          Potential off-target list
 * @param[in]       onTargets           On-target list
 * @param[in]       snpInfoTable        Contains all parts of the fasta ID of SNP sequences
 * @param[in]       seqLookUpTable      List containing [begin, end) of single fasta sequences in combined fasta
 *                                      sequence
 * @param[in]       seqLength           Length of the target sequences that were searched
 */
void filterSnpAlignment(std::vector<unsigned> & validSnpAlignment, std::vector<PotentialOffTarget> & offTargets,
std::map<CharString, PotentialOffTarget> const & onTargets, std::vector<StringSet<CharString> > const & snpInfoTable,
unsigned const & seqLength)
{
    validSnpAlignment.reserve(offTargets.size());

    // Restore original start positions and delete duplicates and on-targets
    for (unsigned i = 0; i < offTargets.size(); i++)
    {
        bool valid = true;

        // Adapt chromosome, position and tag
        StringSet<CharString> fastaID;
        strSplit(fastaID, offTargets[i].chr, EqualsChar<'_'>());
        offTargets[i].chr = fastaID[0];
        offTargets[i].pos = offTargets[i].pos + std::atoi(toCString(fastaID[1]));

        // Check if on-target
        if (comp(offTargets[i], onTargets.at(offTargets[i].target)))
        {
            valid = false;
        }

        getSnpType(offTargets[i].snpType, fastaID, offTargets[i].pos, seqLength);

        // Delete duplicate entries (must be adjacent as entries are sorted)
        // Might happen if off-target lies in REF-only region of 2 variant sequences
        // Can only happen for indels as they have one position too much
        if (i > 0 && comp(offTargets[i], offTargets[i - 1]))
        {
            valid = false;
        }

        if (valid)
        {
            validSnpAlignment.push_back(i);
        }
    }
}

/*!
 * @fn getMismatchPositions
 * @brief Get positions of mismatches in potential off-target based on MD tag in SAM/BAM format
 *
 * @signature void getMismatchPositions(mismatchPos, mismatchString, numMismatches)
 *
 * @param[in,out]   mismatchPos         Positions of mismatches in potential off-target
 * @param[in]       mismatchString      MD tag from SAM/BAM record
 * @param[in]       numMismatches       Maximum number of mismatches allowed during preceeding read mapping
 * @param[in]       strand              Strand of potential off-target
 */
void getMismatchPositions(std::vector<int> & mismatchPos, CharString const & mismatchString,
unsigned const & numMismatches, char & strand)
{
    mismatchPos.reserve(numMismatches); // TODO Change that
    std::istringstream is((std::string) toCString(mismatchString));
    char base;
    unsigned num;
    unsigned pos = 0;
    while(is >> num >> base)
    {
        pos += (num + 1);
        mismatchPos.push_back(pos - 1); // Transform to 0-based
    }
    // If perfect match fill with -1 for comparison reasons with on-target
    // Comparing uninitialized members might result in problems
    if (mismatchPos.empty())
    {
        mismatchPos.push_back(-1);
    }
}

/*!
 * @fn readBamFile
 * @brief Read in a BAM or SAM file and save information of each record in a PotentialOffTarget object
 *
 * @signature void readBamFile(offTargets, bamFilePath, faiIndex, numMismatches)
 *
 * @param[in,out]   offTargets          List of potential off-targets extracted from file
 * @param[in]       bamFilePath         Path to BAM or SAM file
 * @param[in]       faiIndex            The fai index to be searched
 * @param[in]       numMismatches       Maximum number of mismatches allowed during preceeding read mapping
 */
void readBamFile(std::vector<PotentialOffTarget> & offTargets, char const * bamFilePath, FaiIndex & faiIndex,
unsigned const & numMismatches)
{
    BamFileIn bamFileIn;
    if (!open(bamFileIn, bamFilePath))
    {
        throw std::runtime_error("ERROR: Could not open BAM file.");
    }

    CharString tag = "MD";

    try
    {
        // Read header
        BamHeader header;
        readHeader(header, bamFileIn);

        auto & ioContext = context(bamFileIn);

        BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            PotentialOffTarget pot;
            pot.target = record.qName;
            pot.chr = (std::string) toCString(contigNames(ioContext)[record.rID]);
            pot.pos = record.beginPos;

            const int bit5 = 1 << 4;
            if ((record.flag & bit5) != 0)
            {
                pot.strand = '-';
            }
            else
            {
                pot.strand = '+';
            }
            extractSequenceFromIndex(pot.sequence, faiIndex, pot.chr, pot.pos, pot.pos + 23, pot.strand, false);

            pot.snpType = "REF";

            unsigned tagId;
            BamTagsDict tagsDict(record.tags);
            findTagKey(tagId, tagsDict, tag);
            CharString mismatchString;
            extractTagValue(mismatchString, tagsDict, tagId);
            getMismatchPositions(pot.mismatchPos, mismatchString, numMismatches, pot.strand);
            offTargets.push_back(pot);
        }
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return;
    }

}

/*!
 * @fn getSnpInfoTable
 * @brief Construct a list of fastaID parts as chromosome, position and SNPs (and length of corresponding sequence) and
 * set containing all chromosomes with index (later used in sortedIndexAllChr)
 *
 * @signature void getSnpInfoTable(chrMap, snpInfoTable, ids, seqs)
 *
 * @param[in,out]   chrMap              Set containing all chromosomes of SNP regions and the corresponding index
 *                                      (later used in sortedIndexAllChr)
 * @param[in,out]   snpInfoTable        Contains all parts of the fasta ID of SNP sequences and length of corresponding
 *                                      fasta sequence
 * @param[in]       ids                 IDs from SNP genome fasta file
 * @param[in]       seqs                Sequences from SNP genome fasta file
 */
void getSnpInfoTable(std::map<CharString, unsigned> & chrMap, std::vector<StringSet<CharString> > & snpInfoTable,
StringSet<CharString> const & ids, StringSet<Dna5String> const & seqs)
{
    snpInfoTable.resize(length(ids));
    for (unsigned i = 0; i < length(ids); i++)
    {
        StringSet<CharString> fastaID;
        strSplit(fastaID, ids[i], EqualsChar<'_'>());
        reserve(fastaID, length(fastaID) + 1);
        appendValue(fastaID, (CharString) std::to_string(length(seqs[i])));

        chrMap.insert(std::make_pair(fastaID[0], chrMap.size()));
        snpInfoTable[i] = fastaID;
    }
}

/*!
 * @fn readOntargets
 * @brief Read in a BED file and save information of each record in a PotentialOffTarget object
 *
 * @signature void readOntargets(onTargets, offTargetCount, ontargetPath, faiIndex, numMismatches)
 *
 * @param[in,out]   offTargets          Map of on-targets
 * @param[in,out]   offTargetCount      Map of on-targets with count for off-targets
 * @param[in]       ontargetPath        Path to on-target BED6 file
 * @param[in]       faiIndex            The fai index to be searched
 * @param[in]       numMismatches       Maximum number of mismatches allowed during preceeding read mapping
 */
void readOntargets(std::map<CharString, PotentialOffTarget> & onTargets,
std::map<CharString, unsigned> & offTargetCount, char const * ontargetPath, FaiIndex & faiIndex,
unsigned const & numMismatches)
{
    BedFileIn bedIn;

    if (!open(bedIn, ontargetPath))
    {
        throw std::runtime_error("ERROR: Could not open BED file.");
    }

    BedRecord<Bed6> record;
    try
    {
        while (!atEnd(bedIn))
        {
            readRecord(record, bedIn);
            PotentialOffTarget pot;
            pot.target = record.name;
            pot.chr = record.ref;
            pot.pos = record.beginPos;
            pot.strand = record.strand;
            extractSequenceFromIndex(pot.sequence, faiIndex, pot.chr, pot.pos, pot.pos + 23, pot.strand, false);
            pot.mismatchPos = std::vector<int>{-1};
            pot.snpType = "REF";
            onTargets.insert(std::make_pair(record.name, pot));
            offTargetCount.insert(std::make_pair(record.name, 0));
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }
}

}
