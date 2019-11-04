#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

#include "common.h"

using namespace std;
using namespace seqan;

int main(int argc, char *argv[])
{
    typedef StringSet<Dna5String, Owner<ConcatDirect<> > >  TText;
    typedef Index<TText, TIndexConfig>                      TIndex;

    // Argument Parser
    ArgumentParser parser("VARSCOT - Index Creation");
    addDescription(parser, "Application for creating an index (that is needed for the search application). Only supports multi sequence FASTA-files with Dna5 alphabet (A, C, G, T, N). The FASTA file may not contain more than 4 giga bases in total. The index is built using secondary memory. If you get an IO-error, you are probably running out of quota. Try: TMPDIR=/path/to/somewhere/with/enough/quota");

    addOption(parser, ArgParseOption("G", "genome", "Path to the genome (.fa, .fasta, .fastq)", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "genome", "fa fasta fastq");
	setRequired(parser, "genome");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "index");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString genomePath, indexPath;
    getOptionValue(genomePath, parser, "genome");
    getOptionValue(indexPath, parser, "index");

    // Read fasta input file
    StringSet<CharString> ids;
    TText genome;
    SeqFileIn seqFileIn(toCString(genomePath));
    readRecords(ids, genome, seqFileIn);
    clear(ids);

    cout << "Number of sequences: " << length(genome) << endl;

    // Build and save index
    TIndex index(genome);
    indexCreate(index, FibreSALF());
    save(index, toCString(indexPath));

    cout << "Index created successfully" << endl;

    return 0;
}
