#include <iostream>

#include "extract_fasta_ontargets.h"


using namespace seqan;

int main(int argc, char const ** argv)
{

    if (argc != 5)
    {
        std::cerr << "USAGE: extract_fasta_ontargets OUTPUT1.fa OUTPUT2.fa ONTARGETS.bed GENOME.fa\n";
        return 1;
    }

    // Fasta without flanking regions for read mapping
    try
    {
        writeFastaOntargets(argv[1], argv[3], argv[4], false);
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    // Fasta with flanking regions for TUSCAN
    try
    {
        writeFastaOntargets(argv[2], argv[3], argv[4], true);
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}

