#include <iostream>
#include <vector>

#include "process_vcf.h"
#include "overlap_sequences.h"
#include "write_fasta.h"


using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 7)
    {
        std::cerr << "USAGE: vcf_loader FILE.vcf SNPGENOME.fa GENOME.fa SAMPLE SEQLENGTH THREADS\n";
        return 1;
    }

    unsigned sample;
    if (!lexicalCast(sample, argv[4]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[4] << " into an unsigned.\n";
        return 1;
    }

    unsigned seqLength;
    if (!lexicalCast(seqLength, argv[5]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[5] << " into an unsigned.\n";
        return 1;
    }

    unsigned threads;
    if (!lexicalCast(threads, argv[6]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[6] << " into an unsigned.\n";
        return 1;
    }

    // Read in vcf file
    std::cout << "Process records" << std::endl;
    std::vector<std::vector<VariantSequence> > allVariants;
    std::vector<CharString> chrTable;

    try
    {
        processVcfFile(allVariants, chrTable, argv[1], seqLength, sample);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    // Compute overlapping variant ranges
    unsigned chrNumber = chrTable.size();
    std::vector<std::vector<unsigned> > sortedIndex;
    std::vector<std::vector<Pair<unsigned, unsigned> > > allOverlapRegions;
    std::vector<std::vector<unsigned> > indexCenterVariants;
    std::cout << "Compute overlap sequences" << std::endl;
    getVariantOverlapRanges(allOverlapRegions, indexCenterVariants, sortedIndex, chrNumber, allVariants, seqLength, threads);

    // Write fasta file
    std::cout << "Write fasta" << std::endl;

    try
    {
        writeFastaFile(argv[2], argv[3], allVariants, sortedIndex, allOverlapRegions, indexCenterVariants, chrTable, seqLength);
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
