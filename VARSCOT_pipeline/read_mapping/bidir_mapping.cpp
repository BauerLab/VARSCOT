#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>

#include "common.h"

using namespace std;
using namespace seqan;

typedef Pair<uint16_t, uint32_t> TOccType;

struct BamRecord
{
    BamAlignmentRecord r;
    unsigned mismatches;
};

bool isValidPAM(Dna5String const & pam, std::vector<Dna5String> const & validPAM)
{
    for (auto pamF : validPAM)
    {
        if (pam == pamF)
            return true;
    }
    return false;
}

template <typename TIndex, typename TMap>
void searchAndVerify(TIndex & index, CharString const & id, DnaString const & fullRead, DnaString const & partialRead, bool const reverseStrand, bool const firstHalf, TMap & records, unsigned const maxMismatches, std::vector<Dna5String> const & validForwardPAM, std::vector<Dna5String> const & validReversePAM)
{
    auto delegate = [& index, & id, & fullRead, & partialRead, reverseStrand, firstHalf, & records, & maxMismatches, & validForwardPAM, & validReversePAM](auto const & it, DnaString const & /*needle*/, unsigned const /*errors*/
        #ifdef ENABLE_DEBUG_MACRO
            , string str
        #endif
    ) {
        for (auto occ : getOccurrences(it))
        {
            unsigned mismatches = 0;
            auto chromosomeId = occ.i1;
            auto posInChromosome = occ.i2;
            Dna5String mappedRegion;

            auto const & chromosome = indexText(index)[chromosomeId];

            if (firstHalf) // Extend to the right
            {
                // check whether there is enough space left in the genome to the RIGHT, otherwise: continue;
                if (length(chromosome) <= posInChromosome + length(fullRead))
                    continue;
            }
            else // Extend to the left
            {
                // check whether there is enough space left in the genome to the LEFT, otherwise: continue;
                if ((signed) posInChromosome - (signed) (length(fullRead) - length(partialRead)) < 0)
                    continue;

                // make posInChrome point to position where fullRead begins...
                posInChromosome -= length(fullRead) - length(partialRead);
            }

            if (records.find(TOccType(chromosomeId, posInChromosome)) != records.end())
                continue;

            mappedRegion = infixWithLength(chromosome, posInChromosome, length(fullRead));
            //cout << "Mapped Region in Genome: " << mappedRegion << endl;

            // does it end with GG or any other allowed PAM?
            if (!reverseStrand && !isValidPAM(suffix(mappedRegion, length(mappedRegion) - 2), validForwardPAM))
                continue;

            // does it start with CC or any other allowed PAM reverse?
            if (reverseStrand && !isValidPAM(prefix(mappedRegion, 2), validReversePAM))
                continue;

            // at most 8 mismatches?
            for (unsigned i = 0; i < length(fullRead) && mismatches <= maxMismatches; ++i)
            {
                if (ordEqual(chromosome[posInChromosome + i], Dna5('N')))
                    mismatches += maxMismatches + 1; // make alignment "invalid"
                mismatches += !ordEqual(fullRead[i], chromosome[posInChromosome + i]);
            }
            if (mismatches > maxMismatches)
                continue;

            BamRecord record;
            record.r.qName = id;

            // We will get 4 different flag values in the end
            // 0 if forward strand and best match
            // 16 if reverse strand and best match
            // 256 if forward strand and secondary alignment
            // 272 if reverse strand and secondary alignment
            // 272 if reverse strand and secondary alignment
            if (reverseStrand)
                record.r.flag = BAM_FLAG_RC;
            record.r.rID = chromosomeId;
            record.r.beginPos = posInChromosome;
            record.r.mapQ = 255;
            appendValue(record.r.cigar, CigarElement<>('M', 23));
            record.r.rNextId = BamAlignmentRecord::INVALID_REFID;
            record.r.pNext = BamAlignmentRecord::INVALID_POS;
            record.r.tLen = BamAlignmentRecord::INVALID_LEN;
            record.r.seq = fullRead;
            if (reverseStrand)
                reverseComplement(record.r.seq);
            record.r.qual = "IIIIIIIIIIIIIIIIIIIIIII";

            CharString rawTagsText;
            BamTagsDict tagsDict(rawTagsText);

            CharString md;
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), mappedRegion);
            assignSource(row(align, 1), fullRead);
            getMDString(md, row(align, 0), row(align, 1));

            setTagValue(tagsDict, "NM", mismatches);
            setTagValue(tagsDict, "MD", md);
            record.r.tags = host(tagsDict);
            record.mismatches = mismatches;
            records.insert(pair<TOccType, BamRecord>(TOccType(chromosomeId, posInChromosome), record));
        }
    };

    switch (maxMismatches)
    {
        case 0: case 1:
            find<0, 0>(delegate, index, partialRead, HammingDistance());
            break;
        case 2: case 3:
            find<0, 1>(delegate, index, partialRead, HammingDistance());
            break;
        case 4: case 5:
            find<0, 2>(delegate, index, partialRead, HammingDistance());
            break;
        case 6: case 7:
            find<0, 3>(delegate, index, partialRead, HammingDistance());
            break;
        case 8:
            find<0, 4>(delegate, index, partialRead, HammingDistance());
            break;
    }

}

template <typename TIndex, typename TContext>
void searchAndVerifyEntireRead(TIndex & index, CharString const & id, DnaString const & read, bool const reverseStrand, TContext const & bamIOContext,
String<char> & buffer, unsigned const mismatches, std::vector<Dna5String> const & validForwardPAM, std::vector<Dna5String> const & validReversePAM)
{
    map<TOccType, BamRecord> records;

    // erste Hälfte
    DnaString firstHalf(infix(read, 0, length(read)/2));
    searchAndVerify(index, id, read, firstHalf, reverseStrand, true, records, mismatches, validForwardPAM, validReversePAM);

    // zweite Hälfte
    DnaString secondHalf(infix(read, length(read)/2, length(read)));
    searchAndVerify(index, id, read, secondHalf, reverseStrand, false, records, mismatches, validForwardPAM, validReversePAM);

    if (records.empty())
        return;

    auto best_record_pair = records.begin();

    // set secondary flag
    for (auto it = next(records.begin()); it != records.end(); ++it)
    {
        auto & record_pair = *it;

        if (record_pair.second.mismatches >= (*best_record_pair).second.mismatches) // not less mismatches
        {
            record_pair.second.r.flag |= BAM_FLAG_SECONDARY;
            write(buffer, record_pair.second.r, bamIOContext, Sam());
        }
        else
        {
            (*best_record_pair).second.r.flag |= BAM_FLAG_SECONDARY;
            write(buffer, (*best_record_pair).second.r, bamIOContext, Sam());
            best_record_pair = it;
        }
    }

    write(buffer, (*best_record_pair).second.r, bamIOContext, Sam());
}

int main(int argc, char *argv[])
{
    // Argument parser
    ArgumentParser parser("Read mapping");
    addDescription(parser, "Read mapper for CRISPR-Cas9 off-targets with a bidirectional FM index. Only supports Dna4 so far (everything else than ACGT will be converted to A). All reads must have the same length.");

    addOption(parser, ArgParseOption("G", "genome", "Path to genome fasta file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "genome", "fa fasta");
    setRequired(parser, "genome");

    addOption(parser, ArgParseOption("I", "index", "Path to the indexed genome", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");

    addOption(parser, ArgParseOption("R", "reads", "Path to the reads (have to be Dna4)", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "reads", "fa fasta");
    setRequired(parser, "reads");

    addOption(parser, ArgParseOption("M", "mismatches", "Number of allowed mismatches", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "mismatches");

    addOption(parser, ArgParseOption("T", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", "1");

    addOption(parser, ArgParseOption("O", "output", "Path to output SAM file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setValidValues(parser, "output", "sam bam");

    addOption(parser, ArgParseOption("P", "pam", "Additional non-canonical PAM that should be allowed for off-target search besides (N)GG and (N)GA (default).", seqan::ArgParseArgument::STRING, "STR"));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    unsigned mismatches, threads;
    CharString genomePath, indexPath, readsPath, outputPath;
    Dna5String additionalPAM;
    getOptionValue(genomePath, parser, "genome");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(readsPath, parser, "reads");
    getOptionValue(mismatches, parser, "mismatches");
    getOptionValue(threads, parser, "threads");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(additionalPAM, parser, "pam");

	if (mismatches < 0 || mismatches > 8)
	{
		cerr << "Error: Maximum number of mismatches must lie between 0 and 8." << endl;
		return 1;
	}

    std::vector<Dna5String> validForwardPAM{"GG", "GA"};
    std::vector<Dna5String> validReversePAM{"CC", "TC"};
    if (additionalPAM != "")
    {
        validForwardPAM.push_back(additionalPAM);
        reverseComplement(additionalPAM);
        validReversePAM.push_back(additionalPAM);
    }

    // Index configuration
    typedef StringSet<Dna5String, Owner<ConcatDirect<> > > TText;
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
    typedef Index<TText, BidirectionalIndex<FMIndex<void, TMyFastConfig> > > TIndex;

    TIndex index;
    StringSet<CharString> ids;
    StringSet<DnaString> reads;
    vector<String<char> > output_buffer;
    output_buffer.resize(threads);
    StringSet<CharString> contigNameStore;

    // Open reads
    // TODO: throw errror when Dna5String
    SeqFileIn seqFileIn(toCString(readsPath));
    readRecords(ids, reads, seqFileIn);
    cout << "Reads loaded (total: " << length(reads) << ")." << endl;

    // Open index
    open(index, toCString(indexPath));
    cout << "Index loaded." << endl;

    // Read ids from genome fasta: TODO improve
    SeqFileIn genomeFileIn(toCString(genomePath));
    while (!atEnd(genomeFileIn))
    {
        CharString id;
        Dna5String genome;
        readRecord(id, genome, genomeFileIn);
        appendValue(contigNameStore, id);
        // throw away seq
    }

    NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);

    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (unsigned i = 0; i < length(reads); ++i)
    {
        DnaString read(infix(reads[i], 0, length(reads[i])));
        // Search & verify reads
        searchAndVerifyEntireRead(index, ids[i], read, false, bamIOContext, output_buffer[omp_get_thread_num()], mismatches, validForwardPAM, validReversePAM);
        // Search & verify reverse strand!
        reverseComplement(read);
        searchAndVerifyEntireRead(index, ids[i], read, true, bamIOContext, output_buffer[omp_get_thread_num()], mismatches, validForwardPAM, validReversePAM);
    }

    // TODO: SeqAn SAM/BAM IO?
    ofstream out;
    out.open(toCString(outputPath));

    if (!out.is_open())
    {
        cerr << "ERROR: Could not open output path." << endl;
		return 1;
    }

    for (auto i : output_buffer)
        out << i;
    out.close();

    return 0;
}
