// ============================================================================
// Extract on-target sequences from fai index and write them to a fasta file
// ============================================================================

#include <iostream>

#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#pragma once

namespace seqan
{

/*!
 * @fn extractSequenceFromIndex
 * @brief Extracts a sequence from an fai index based on chromosome, start, end and strand (possible reverse
 * complement). Extracts flanking regions (total of 30 bp) if specified.
 *
 * @signature void extractSequenceFromIndex(sequence, faiIndex, id, start, end)
 *
 * @param[in,out]   sequence        Extracted sequence
 * @param[in]       faiIndex        The fai index to be searched
 * @param[in]       id              The fasta ID were the sequence is located
 * @param[in]       start           Start position of the sequence (0-based)
 * @param[in]       end             End position of the sequence (not included, half open intervals)
 * @param[in]       strand          Strand of sequence (if reverse complement needed)
 * @param[in]       flanking        Bool defining whether flanking regions should be extracted
 *
 * @throw Exception if sequence cannot be read in
 */
void extractSequenceFromIndex(Dna5String & sequence, FaiIndex & faiIndex, CharString & id, unsigned start, unsigned end,
char & strand, bool flanking)
{
    // Translate sequence name to index
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, id))
    {
        throw std::out_of_range("ERROR: Index out of range.");
    }

    // Check flanking/no flanking option and adjust start and end if needed
    if (flanking && strand == '+')
    {
        start -= 4;
        end += 3;
    }
    else if (flanking && strand == '-')
    {
        start -= 3;
        end += 4;
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
        if (strand == '-')
        {
            reverseComplement(sequence);
        }
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return;
    }
}

/*!
 * @fn writeFastaFile
 * @brief Writes variant sequences which are created on the fly to a fasta file
 *
 * @signature void writeFastaFile(outputPath, indexPath, allVariants, sortedIndexAllChr, allOverlapRegions, chrTable,
 *                                seqLength)
 *
 * @param[in]       outputPath      Output path for fasta file
 * @param[in]       inputPath       Path for input bed6 file
 * @param[in]       indexPath       Path for the genome index should be built of or index exists in same directory
 * @param[in]       flanking        Bool defining whether flanking regions should be extracted
 *
 * @throw Exception if bed records cannot be read in
*/
void writeFastaOntargets(char const * outputPath, char const * inputPath, char const * indexPath, bool flanking)
{
    SeqFileOut seqOut;

    if (!open(seqOut, outputPath))
    {
        throw std::runtime_error("ERROR: Could not open output file.");
    }

    FaiIndex faiIndex;
    if (!open(faiIndex, indexPath))
    {
        if (!build(faiIndex, indexPath))
        {
            throw std::runtime_error("ERROR: Index could not be loaded or built.");
        }
        if (!save(faiIndex))
        {
            throw std::runtime_error("ERROR: Index could not be written do disk.");
        }
    }

    BedFileIn bedIn;

    if (!open(bedIn, inputPath))
    {
        throw std::runtime_error("ERROR: Could not open BED file.");
    }

    // Read the bed file
    BedRecord<Bed6> record;
    try
    {
        while (!atEnd(bedIn))
        {
            readRecord(record, bedIn);
            Dna5String sequence;
            extractSequenceFromIndex(sequence, faiIndex, record.ref, record.beginPos, record.endPos, record.strand,
                                     flanking);
            writeRecord(seqOut, record.name, sequence);
        }
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return;
    }
}
}
