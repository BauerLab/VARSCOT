#include <iostream>

#include "merge_output_bam.h"


using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 13)
    {
        std::cerr << "USAGE: bam_merger RESULT_MERGED.txt FEATURE_MATRIX.txt RESULT_REF.bam RESULT_SNP.bam ONTARGETS.bed GENOME.fa VARIANT_GENOME.fa TUSCAN_REGRESSION.txt NUMMISMATCHES SEQLENGTH THREADS MIT\n";
        return 1;
    }

    unsigned numMismatches;
    if (!lexicalCast(numMismatches, argv[9]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[9] << " into an unsigned.\n";
        return 1;
    }

    unsigned seqLength;
    if (!lexicalCast(seqLength, argv[10]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[10] << " into an unsigned.\n";
        return 1;
    }

    unsigned threads;
    if (!lexicalCast(threads, argv[11]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[11] << " into an unsigned.\n";
        return 1;
    }

    unsigned mit;
    if (!lexicalCast(mit, argv[12]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[12] << " into an unsigned.\n";
        return 1;
    }

    try
    {
        if (mit == 0)
        {
            mergeResults(argv[1], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], numMismatches, seqLength, threads);
        }
        else
        {
            mergeResults(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], numMismatches, seqLength, threads);
        }
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
