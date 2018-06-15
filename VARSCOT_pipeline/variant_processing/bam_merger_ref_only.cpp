#include <iostream>

#include "merge_output_bam.h"


using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 10)
    {
        std::cerr << "USAGE: bam_merger RESULT_MERGED.txt FEATURE_MATRIX.txt RESULT_REF.bam ONTARGETS.bed GENOME.fa TUSCAN_REGRESSION.txt NUMMISMATCHES SEQLENGTH MIT\n";
        return 1;
    }

    unsigned numMismatches;
    if (!lexicalCast(numMismatches, argv[7]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[7] << " into an unsigned.\n";
        return 1;
    }

    unsigned seqLength;
    if (!lexicalCast(seqLength, argv[8]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[8] << " into an unsigned.\n";
        return 1;
    }

    unsigned mit;
    if (!lexicalCast(mit, argv[9]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[9] << " into an unsigned.\n";
        return 1;
    }

    try
    {
        if (mit == 0)
        {
            processRefOnly(argv[1], argv[3], argv[4], argv[5], argv[6], numMismatches, seqLength);
        }
        else
        {
            processRefOnly(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], numMismatches, seqLength);
        }
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
