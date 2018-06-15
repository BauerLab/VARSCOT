// ============================================================================
// Read in a VCF file and process each record
// ============================================================================
#pragma once

#include <iostream>
#include <vector>

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>
#include <sstream>

namespace seqan
{

/*!
 * @struct VariantSequence
 * @brief Stores information needed to define a variant sequence
 *
 * @signature struct VariantSequence;
 *
 * The VariantSequence object stores information about a single variant (substitution, insertion or deletion).
 * This information includes chromosome (index of reading order of chromosomes), position of the first base of the
 * variant, end of the sequence that will be generated, reference, alternative and variant type.
 * Additionally the genotype information is saved. If the variant is phased, information about the allele the variant
 * lies on is stored (0 for mother, 1 for father, 2 for both and -1 for unphased).
 * Note: Positions are 0-based. As half open intervals are used end position is not included in the later sequence.
 */
struct VariantSequence
{
    Dna5String ref;
    Dna5String alt;
    unsigned chr;
    unsigned pos;
    unsigned start;
    unsigned end;
    unsigned variantType; // 0 for substitution, 1 for insertion, 2 for deletion
    int allele; // 0 for first (mother), 1 for second (father), 2 for both and -1 for unphased
};

/*!
 * @fn processRecord
 * @brief Processes a single record of a VCF file and stores every included variant in a VariantSequence object
 *
 * @signature void processRecord(variants, record, seqLength)
 *
 * @param[in,out]   variants        Vector containing all variant objects from the VCF record
 * @param[in]       record          VCF record that is processed by the function
 * @param[in]       seqLength       Length of the target sequence that should be searched later (normally 23 bp)
 * @param[in]       sampleIndex     Index of the sample that shoud be processed (0-based)
 */
void processRecord(std::vector<VariantSequence> & variants, VcfRecord & record, unsigned const & seqLength, unsigned const & sampleIndex)
{
    VariantSequence vs;
    vs.chr = record.rID;
    vs.pos = record.beginPos;
    vs.ref = record.ref;

    if (sampleIndex >= length(record.genotypeInfos))
    {
        throw std::out_of_range("ERROR: Sample index out of range.");
    }

    StringSet<CharString> individualGenotypeInfosEntries;
    strSplit(individualGenotypeInfosEntries, record.genotypeInfos[sampleIndex], EqualsChar<':'>());

    // Check which position contains the genotype information
    StringSet<CharString> formatEntries;
    strSplit(formatEntries, record.format, EqualsChar<':'>());

    unsigned positionGT;
    CharString tagGT = "GT";

    for (unsigned i = 0; i < length(formatEntries); i++)
    {
        if (formatEntries[i] == tagGT)
        {
            positionGT = i;
            break;
        }
    }
    int firstAllele = -1;
    int secondAllele = -1;
    char sep;
    bool phased = true;

    // Split alt column into single alternative sequences
    StringSet<CharString> altAlleles;
    strSplit(altAlleles, record.alt, EqualsChar<','>());

    std::istringstream is((std::string) toCString(individualGenotypeInfosEntries[positionGT]));

    if (is >> firstAllele && firstAllele <= length(altAlleles))
    {
        if (is >> sep >> secondAllele && secondAllele <= length(altAlleles))
        {
            if (sep == '/')
            {
                phased = false;
            }
        }
        else
        {
            // This case exists only for the Y chromosome when there exist only one allele
            secondAllele = firstAllele;
        }
    }
    else
    {
        return;
    }

    // Handle all possible cases - phased and unphased
    if (firstAllele == 0u && secondAllele == 0u)
    {
        // 0|0  only reference
        return;
    }
    else if (firstAllele > 0u && secondAllele > 0u && firstAllele != secondAllele)
    {
        // Only case where there exist 2 variants
        // Check if missing alternative values exist
        if (altAlleles[firstAllele - 1] != '.' && altAlleles[secondAllele - 1] != '.')
        {
            variants.resize(2);

            vs.allele = 0;
            vs.alt = altAlleles[firstAllele - 1];
            variants[0] = vs;

            vs.allele = 1;
            vs.alt = altAlleles[secondAllele - 1];
            variants[1] = vs;
        }
        else if (altAlleles[firstAllele - 1] != '.')
        {
            variants.resize(1);

            vs.allele = 0;
            vs.alt = altAlleles[firstAllele - 1];
            variants[0] = vs;
        }
        else if (altAlleles[secondAllele - 1] != '.')
        {
            variants.resize(1);

            vs.allele = 1;
            vs.alt = altAlleles[secondAllele - 1];
            variants[1] = vs;
        }
        else
        {
            return;
        }
    }
    else
    {
        if (altAlleles[0] == '.')
        {
            return;
        }

        variants.resize(1);

        if (firstAllele == 0u)
        {
            // 0|1, 0|2, ...  one reference, one alternative
            vs.allele = 1;
            vs.alt = altAlleles[secondAllele - 1];
        }
        else if (secondAllele == 0u)
        {
            // 1|0, 2|0, ...  one alternative, one reference
            vs.allele = 0;
            vs.alt = altAlleles[firstAllele - 1];
        }
        else
        {
            // 1|1, 1/1, 2|2, 2/2, ... both same alternative
            vs.allele = 2;
            vs.alt = altAlleles[firstAllele - 1]; // no difference between first or second
        }
        variants[0] = vs;
    }

    // Create object for every alternative
    for (unsigned i = 0; i < variants.size(); i++)
    {
        if (phased == false && firstAllele != secondAllele)
        {
            variants[i].allele = -1;
        }
        // Window needs to be calculated from largest entry (either ref or alt)
        if (length(variants[i].ref) > length(variants[i].alt))
        {
            variants[i].variantType = 2; // deletion
        }
        else if (length(variants[i].ref) == length(variants[i].alt))
        {
            variants[i].variantType = 0; // substitution
        }
        else
        {
            variants[i].variantType = 1; // insertion
        }
    }
}

/*!
 * @fn processVcfFile
 * @brief Processes a single record of a VCF file and stores every included variant in a VariantSequence object
 *
 * @signature void processVcfFile(allVariants, chrTable, inputPath)
 *
 * @param[in,out]   allVariants     Vector of vectors containing all variants from the complete VCF file, variants from
 *                                  the same record are stored in the same vector
 * @param[in,out]   chrTable        Look up table containing the chromosome names at the corresponding index
 * @param[in]       inputPath       Path to the VCF file
 * @param[in]       seqLength       Length of the target sequence that should be searched later (normally 23 bp)
 * @param[in]       sampleIndex     Index of the sample that shoud be processed (0-based)
 *
 * @throw Exception if records cannot be read in
 */
void processVcfFile(std::vector<std::vector<VariantSequence> > & allVariants, std::vector<CharString> & chrTable,
char const * inputPath, unsigned const & seqLength, unsigned const & sampleIndex)
{
    // Open input file
    VcfFileIn vcfIn;

    if (!open(vcfIn, inputPath))
    {
        throw std::runtime_error("ERROR: Could not open VCF file.");
    }

    VcfHeader header;
    VcfRecord record;

    try
    {
        // Copy over header
        readHeader(header, vcfIn);

        // Copy the file record by record
        while (!atEnd(vcfIn))
        {
            readRecord(record, vcfIn);
            std::vector<VariantSequence> vs;
            processRecord(vs, record, seqLength, sampleIndex);
            if (!vs.empty())
                allVariants.push_back(vs);
        }
        // Read out the io context and save in look up table
        // This is needed because chromosome names are stored as index of reading order
        // Example: First VCF record has chr20, so 0 is saved instead of chr20
        // Look up table to benefit of this index and later get correct chromosome name
        auto & ioContext = context(vcfIn);
        chrTable.resize(length(contigNames(ioContext)));

        for (unsigned i = 0; i < length(contigNames(ioContext)); i++)
            chrTable[i] = "chr" + (std::string) toCString(contigNames(ioContext)[i]);
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return;
    }
}

}
