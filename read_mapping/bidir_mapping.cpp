#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>

#include "helper.h"
#include "searchSchemes.h"

using namespace std;
using namespace seqan;

namespace seqan {

template <typename TChar, typename TOwner>
struct SAValue<StringSet<String<TChar>, TOwner > >
{
    typedef Pair<uint32_t, uint32_t> Type;
};

template <typename TChar, typename TOwner>
struct SAValue<String<TChar, TOwner > >
{
    typedef uint32_t Type;
};

};

typedef Pair<uint16_t, uint32_t> TOccType;

struct BamRecord
{
    BamAlignmentRecord r;
    unsigned mismatches;
};

template <typename TIndex, typename TMap>
void searchAndVerify(TIndex & index, CharString const & id, DnaString const & fullRead, DnaString const & partialRead, bool const reverseStrand, bool const firstHalf,
      TMap & records, unsigned const maxMismatches)
{
    auto scheme = schemes[maxMismatches/2]; // <= 3 errors je Hälfte
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;
    TIter it(index);

    auto delegate = [& index, & id, & fullRead, & partialRead, reverseStrand, firstHalf, & records, & maxMismatches](auto const & it, unsigned const errors
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

            // does it end with GG?
            if (!reverseStrand && suffix(mappedRegion, length(mappedRegion) - 2) != Dna5String("GG") && suffix(mappedRegion, length(mappedRegion) - 2) != Dna5String("GA"))
                continue;

            // does it start with CC?
            if (reverseStrand && prefix(mappedRegion, 2) != Dna5String("CC") && prefix(mappedRegion, 2) != Dna5String("TC"))
                continue;

            // at most 8 mismatches?
            // Sara: Wenn ich eine Haelfte mit 3 und eine Haelfte mit 5 Fehlern suche, dann kriege ich nicht die Ergebnisse die 4 msimatches in jedem Teil haben oder?
            for (unsigned i = 0; i < length(fullRead) && mismatches <= maxMismatches; ++i)
            {
                if (ordEqual(chromosome[posInChromosome + i], Dna5('N')))
                    mismatches += (maxMismatches + 1); // make alignment "invalid"
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
            {
                record.r.flag = BAM_FLAG_RC;
            }
            record.r.rID = chromosomeId;
            record.r.beginPos = posInChromosome;
            record.r.mapQ = 255;
            appendValue(record.r.cigar, CigarElement<>('M', 23));
            record.r.rNextId = BamAlignmentRecord::INVALID_REFID;
            record.r.pNext = BamAlignmentRecord::INVALID_POS;
            record.r.tLen = BamAlignmentRecord::INVALID_LEN;
            record.r.seq = fullRead;
            if (reverseStrand)
            {
                reverseComplement(record.r.seq);
            }
            record.r.qual = "IIIIIIIIIIIIIIIIIIIIIII";

            // Sara: Added NM and MD tag
            // Sara: For MD string we need Gap Objects - This is an ugly fix by converting them to align - better idea?
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

    computeBlocklength<typename Value<DnaString>::Type>(scheme, length(partialRead), false /* optimalBlocklength */);
    for (Search const & s : scheme)
    {
        goRoot(it);
        search(delegate, it, partialRead, s, false /*indels*/, false /*debug*/);
    }
}

template <typename TIndex, typename TContext>
void searchAndVerifyEntireRead(TIndex & index, CharString const & id, DnaString const & read, bool const reverseStrand,
      TContext const & bamIOContext, String<char> & buffer, unsigned const mismatches)
{
    map<TOccType, BamRecord> records;

    // erste Hälfte
    DnaString firstHalf(infix(read, 0, length(read)/2));
    searchAndVerify(index, id, read, firstHalf, reverseStrand, true, records, mismatches);

    // zweite Hälfte
    DnaString secondHalf(infix(read, length(read)/2, length(read)));
    searchAndVerify(index, id, read, secondHalf, reverseStrand, false, records, mismatches);

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
    ArgumentParser parser("SearchSchemes - Benchmarking");
    addDescription(parser, "App for benchmarking the running time of search schemes. Only supports Dna4 so far (everything else than ACGT will be converted to A). All reads must have the same length.");

    addOption(parser, ArgParseOption("G", "genome", "Path to genome fasta file", ArgParseArgument::INPUT_FILE, "IN"));
  	setValidValues(parser, "genome", "fa fasta");
    setRequired(parser, "genome");

    addOption(parser, ArgParseOption("I", "index", "Path to the indexed genome", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "index");

    addOption(parser, ArgParseOption("R", "reads", "Path to the reads", ArgParseArgument::INPUT_FILE, "IN"));
  	setValidValues(parser, "reads", "fa fasta");
  	setRequired(parser, "reads");

    addOption(parser, ArgParseOption("M", "mismatches", "Number of allowed mismatches", ArgParseArgument::INTEGER, "INT"));
  	setRequired(parser, "mismatches");

  	addOption(parser, ArgParseOption("T", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
  	setRequired(parser, "threads"); // TODO: Necessary? How to set a default value?

    addOption(parser, ArgParseOption("O", "output", "Path to output SAM file", ArgParseArgument::OUTPUT_FILE, "OUT"));
  	setValidValues(parser, "output", "sam bam");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    unsigned mismatches, threads;
    CharString genomePath, indexPath, readsPath, outputPath;
    getOptionValue(genomePath, parser, "genome");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(readsPath, parser, "reads");
    getOptionValue(mismatches, parser, "mismatches");
    getOptionValue(threads, parser, "threads");
    getOptionValue(outputPath, parser, "output");

	if (mismatches < 0 || mismatches > 8)
	{
		std::cerr << "Error: Maximum number of mismatches must lie between 0 and 8." << std::endl;
		return 1;
	}

    // Index configuration
    typedef StringSet<Dna5String, Owner<ConcatDirect<> > > TText;
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
    typedef Index<TText, BidirectionalIndex<FMIndex<void, TMyFastConfig> > > TIndex;

    TIndex index;
    StringSet<CharString> ids;
    StringSet<DnaString> reads;
    std::vector<String<char> > output_buffer;
    output_buffer.resize(threads);
    StringSet<CharString> contigNameStore;

    // Open reads
    // TODO: throw errror when Dna5String
    SeqFileIn seqFileIn(toCString(readsPath));
    readRecords(ids, reads, seqFileIn);
    // cout << "Reads loaded (total: " << length(reads) << ")." << endl;

    // Open index
    open(index, toCString(indexPath));
    // cout << "Index loaded." << endl;

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
	// cout << output_buffer.size() << endl;

	omp_set_num_threads(threads);
	#pragma omp parallel for schedule(static)
    for (unsigned i = 0; i < length(reads); ++i)
    {
		// printf("thread id: %i\n", omp_get_thread_num());
        DnaString read(infix(reads[i], 0, length(reads[i])));
        // Search & verify reads
        searchAndVerifyEntireRead(index, ids[i], read, false, bamIOContext, output_buffer[omp_get_thread_num()], mismatches);
        // Search & verify reverse strand!
        reverseComplement(read);
        searchAndVerifyEntireRead(index, ids[i], read, true, bamIOContext, output_buffer[omp_get_thread_num()], mismatches);
    }

    // TODO: SeqAn SAM/BAM IO?
    std::ofstream out;
    out.open(toCString(outputPath));

    if (out.is_open())
    {
        for (auto i : output_buffer)
            out << i;
        out.close();
    }
    else
    {
        std::cerr << "ERROR: Could not open output path." << std::endl;
		return 1;
    }

    return 0;
}
