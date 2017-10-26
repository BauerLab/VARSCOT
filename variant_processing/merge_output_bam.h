// ============================================================================
// Merge bam output from mapping against reference and SNP genome
// ============================================================================

#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <seqan/parallel.h>
#include <seqan/seq_io.h>

#include "filter_output_bam.h"
#include "feature_matrix.h"
#include "mit_score.h"

namespace seqan
{
/*!
 * @fn mergeResults
 * @brief Merges the SAM/BAM output files of reference and variant genome mapping and produces output with MIT score
 *
 * @signature void mergeResults(outputPath, featureMatrixPath, refBamPath, snpBamPath, ontargetFastaPath,
 *                              indexPathRef, indexPathSnp, tuscanPath, numMismatches, seqLength, threads)
 *
 * @param[in]       outputPath          Path of the merged output file
 * @param[in]       refBamPath          Path to reference genome BAM/SAM file
 * @param[in]       snpFastaPath        Path to SNP genome
 * @param[in]       ontargetPath        Path to ontarget fasta file
 * @param[in]       indexPathRef        Path for the genome, index should be built of or index exists in same directory
 * @param[in]       indexPathSnp        Path for the variant genome, index should be built of or index exists in
 *                                      same directory
 * @param[in]       tuscanPath          Path to on-target activity (TUSCAN regression) file
 * @param[in]       numMismatches       Maximum number of mismatches allowed during preceeding read mapping
 * @param[in]       seqLength           Length of the target sequences that were searched
 * @param[in]       threads             Number of threads used for parallelization
 *
 * This function merges the two output files from mapping agains reference and SNP genome. In addition, it is secured
 * that no duplicates, matches to N region in combined SNP genome, matches to the reference in SNP regions or
 * off-targets with PAM mismatches are included into the final output file.
 * The final output contains for each potential off-target the bed6 records where name is the on-target name, the
 * sequence, the positions of mismatches (if existing) and specific tags if variants are included in the sequence.
 */
void mergeResults(char const * outputPath, char const * refBamPath, char const * snpBamPath, char const * ontargetPath,
char const * indexPathRef, char const * indexPathSnp, char const * tuscanPath, unsigned const & numMismatches,
unsigned const & seqLength, unsigned threads)
{
    // Read in SNP genome
    SeqFileIn fastaFileIn;

    if (!open(fastaFileIn, indexPathSnp))
    {
        throw std::runtime_error("ERROR: Could not open variant genome FASTA file.");
    }

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    try
    {
        readRecords(ids, seqs, fastaFileIn);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
    }

    // Open fai index reference
    FaiIndex faiIndexRef;
    if (!open(faiIndexRef, indexPathRef))
    {
        if (!build(faiIndexRef, indexPathRef))
        {
            throw std::runtime_error("ERROR: Reference index could not be loaded or built.");
        }
        if (!save(faiIndexRef))
        {
            throw std::runtime_error("ERROR: Reference index could not be written do disk.");
        }
    }

    // Open fai index SNP genome
    FaiIndex faiIndexSnp;
    if (!open(faiIndexSnp, indexPathSnp))
    {
        if (!build(faiIndexSnp, indexPathSnp))
        {
            throw std::runtime_error("ERROR: Variant index could not be loaded or built.");
        }
        if (!save(faiIndexSnp))
        {
            throw std::runtime_error("ERROR: Variant index could not be written do disk.");
        }
    }


    // Split fasta IDs for information about SNP sequences
    std::map<CharString, unsigned> chrMap;
    std::vector<StringSet<CharString> > snpInfoTable;
    getSnpInfoTable(chrMap, snpInfoTable, ids, seqs);

    // Read in on-targets and prepare counting number of off-targets per on-target
    std::map<CharString, PotentialOffTarget> onTargets;
    std::map<CharString, unsigned> offTargetCount;

    readOntargets(onTargets, offTargetCount, ontargetPath, faiIndexRef, numMismatches);

    // Prepare for reading and filtering of off-targets
    std::vector<PotentialOffTarget> refOffTargets;
    std::vector<PotentialOffTarget> snpOffTargets;
    std::vector<std::vector<unsigned> > sortedIndexAllChr;

    std::vector<unsigned> validRefAlignment;
    std::vector<unsigned> validSnpAlignment;

    // Process reference and variant genome mapping
    std::cout << "Process reference off-targets" << std::endl;
    readBamFile(refOffTargets, refBamPath, faiIndexRef, numMismatches);
    sortSnpRegionsByChr(sortedIndexAllChr, chrMap, snpInfoTable, threads);
    filterRefAlignment(validRefAlignment, sortedIndexAllChr, chrMap, snpInfoTable, refOffTargets, onTargets, seqLength, threads);

    std::cout << "Process variant off-targets" << std::endl;
    readBamFile(snpOffTargets, snpBamPath, faiIndexSnp, numMismatches);
    filterSnpAlignment(validSnpAlignment, snpOffTargets, onTargets, snpInfoTable, seqLength);

    // Read TUSCAN activity file in
    std::map<CharString, double> onTargetActivity;
    readTuscanResult(onTargetActivity, tuscanPath);

    // Write combined results to a file
    std::ofstream mergedResults;
    mergedResults.open(outputPath);

    if (mergedResults.is_open())
    {
        // Write column names output
        mergedResults << "#Chr\t" << "Start\t" << "End\t" << "Targetsite\t" << "Score\t" << "Strand\t" << "Sequence\t" << "Mismatch_Number\t" << "Mismatch_Positions\t" << "Variants\n";

        // Write reference results
        std::vector<int> perfectMatch{-1};
        double currOntargetActivity;
        std::string offTargetName;

        for (unsigned i = 0; i < validRefAlignment.size(); i++)
        {
            // Count off-targets and get current off-target count for name
            offTargetCount.at(refOffTargets[validRefAlignment[i]].target)++;
            offTargetName = (std::string) toCString(refOffTargets[validRefAlignment[i]].target) + '_' + std::to_string(offTargetCount.at(refOffTargets[validRefAlignment[i]].target));

            mergedResults << refOffTargets[validRefAlignment[i]].chr << "\t";
            mergedResults << refOffTargets[validRefAlignment[i]].pos << "\t";
            mergedResults << (refOffTargets[validRefAlignment[i]].pos + 23) << "\t";
            mergedResults << offTargetName << "\t";
            mergedResults << calcMitScore(refOffTargets[validRefAlignment[i]].mismatchPos) << "\t";
            mergedResults << refOffTargets[validRefAlignment[i]].strand << "\t";
            mergedResults << refOffTargets[validRefAlignment[i]].sequence << "\t";

            if (refOffTargets[validRefAlignment[i]].mismatchPos == perfectMatch)
            {
                mergedResults << 0 << "\t\t";
            }
            else
            {
                mergedResults << refOffTargets[validRefAlignment[i]].mismatchPos.size() << "\t";
                for (unsigned j = 0; j < (refOffTargets[validRefAlignment[i]].mismatchPos.size() - 1); j++)
                {
                    mergedResults << refOffTargets[validRefAlignment[i]].mismatchPos[j] << ',';
                }
                mergedResults << refOffTargets[validRefAlignment[i]].mismatchPos[refOffTargets[validRefAlignment[i]].mismatchPos.size() - 1] << "\t";
            }
            mergedResults << refOffTargets[validRefAlignment[i]].snpType << "\n";
        }

        // Write variant results
        for (unsigned i = 0; i < validSnpAlignment.size(); i++)
        {
            // Count off-targets and get current off-target count for name
            offTargetCount.at(snpOffTargets[validSnpAlignment[i]].target)++;
            offTargetName = (std::string) toCString(snpOffTargets[validSnpAlignment[i]].target) + '_' + std::to_string(offTargetCount.at(snpOffTargets[validSnpAlignment[i]].target));

            mergedResults << snpOffTargets[validSnpAlignment[i]].chr << "\t";
            mergedResults << snpOffTargets[validSnpAlignment[i]].pos << "\t";
            mergedResults << (snpOffTargets[validSnpAlignment[i]].pos + 23) << "\t";
            mergedResults << offTargetName << "\t";
            mergedResults << calcMitScore(snpOffTargets[validSnpAlignment[i]].mismatchPos) << "\t";
            mergedResults << snpOffTargets[validSnpAlignment[i]].strand << "\t";
            mergedResults << snpOffTargets[validSnpAlignment[i]].sequence << "\t";

            if (snpOffTargets[validSnpAlignment[i]].mismatchPos == perfectMatch)
            {
                mergedResults << 0 << "\t\t";
            }
            else
            {
                mergedResults << snpOffTargets[validSnpAlignment[i]].mismatchPos.size() << "\t";
                for (unsigned j = 0; j < (snpOffTargets[validSnpAlignment[i]].mismatchPos.size() - 1); j++)
                {
                    mergedResults << snpOffTargets[validSnpAlignment[i]].mismatchPos[j] << ',';
                }
                mergedResults << snpOffTargets[validSnpAlignment[i]].mismatchPos[snpOffTargets[validSnpAlignment[i]].mismatchPos.size() - 1] << "\t";
            }

            mergedResults << snpOffTargets[validSnpAlignment[i]].snpType << "\n";
        }

        mergedResults.close();
        std::cout << "Merging output files finished" << std::endl;
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open output file.");
    }
}


/*!
 * @fn mergeResults
 * @brief Merges the SAM/BAM output files of reference and variant genome mapping and writes feature matrix as second output
 *
 * @signature void mergeResults(outputPath, featureMatrixPath, refBamPath, snpBamPath, ontargetFastaPath,
 *                              indexPathRef, indexPathSnp, tuscanPath, numMismatches, seqLength, threads)
 *
 * @param[in]       outputPath          Path of the merged output file
 * @param[in]       featureMatrixPath   Path of the merged output file
 * @param[in]       refBamPath          Path to reference genome BAM/SAM file
 * @param[in]       snpFastaPath        Path to SNP genome
 * @param[in]       ontargetPath        Path to ontarget fasta file
 * @param[in]       indexPathRef        Path for the genome, index should be built of or index exists in same directory
 * @param[in]       indexPathSnp        Path for the variant genome, index should be built of or index exists in
 *                                      same directory
 * @param[in]       tuscanPath          Path to on-target activity (TUSCAN regression) file
 * @param[in]       numMismatches       Maximum number of mismatches allowed during preceeding read mapping
 * @param[in]       seqLength           Length of the target sequences that were searched
 * @param[in]       threads             Number of threads used for parallelization
 *
 * This function merges the two output files from mapping agains reference and SNP genome. In addition, it is secured
 * that no duplicates, matches to N region in combined SNP genome, matches to the reference in SNP regions or
 * off-targets with PAM mismatches are included into the final output file.
 * The final output contains for each potential off-target the bed6 records where name is the on-target name, the
 * sequence, the positions of mismatches (if existing) and specific tags if variants are included in the sequence.
 */
void mergeResults(char const * outputPath, char const * featureMatrixPath, char const * refBamPath,
char const * snpBamPath, char const * ontargetPath, char const * indexPathRef,
char const * indexPathSnp, char const * tuscanPath, unsigned const & numMismatches, unsigned const & seqLength,
unsigned threads)
{
    // Read in SNP genome
    SeqFileIn fastaFileIn;

    if (!open(fastaFileIn, indexPathSnp))
    {
        throw std::runtime_error("ERROR: Could not open variant genome FASTA file.");
    }

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    try
    {
        readRecords(ids, seqs, fastaFileIn);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
    }

    // Open fai index reference
    FaiIndex faiIndexRef;
    if (!open(faiIndexRef, indexPathRef))
    {
        if (!build(faiIndexRef, indexPathRef))
        {
            throw std::runtime_error("ERROR: Reference index could not be loaded or built.");
        }
        if (!save(faiIndexRef))
        {
            throw std::runtime_error("ERROR: Reference index could not be written do disk.");
        }
    }

    // Open fai index SNP genome
    FaiIndex faiIndexSnp;
    if (!open(faiIndexSnp, indexPathSnp))
    {
        if (!build(faiIndexSnp, indexPathSnp))
        {
            throw std::runtime_error("ERROR: Variant index could not be loaded or built.");
        }
        if (!save(faiIndexSnp))
        {
            throw std::runtime_error("ERROR: Variant index could not be written do disk.");
        }
    }


    // Split fasta IDs for information about SNP sequences
    std::map<CharString, unsigned> chrMap;
    std::vector<StringSet<CharString> > snpInfoTable;
    getSnpInfoTable(chrMap, snpInfoTable, ids, seqs);

    // Read in on-targets and prepare counting number of off-targets per on-target
    std::map<CharString, PotentialOffTarget> onTargets;
    std::map<CharString, unsigned> offTargetCount;

    readOntargets(onTargets, offTargetCount, ontargetPath, faiIndexRef, numMismatches);

    // Prepare for reading and filtering of off-targets
    std::vector<PotentialOffTarget> refOffTargets;
    std::vector<PotentialOffTarget> snpOffTargets;
    std::vector<std::vector<unsigned> > sortedIndexAllChr;

    std::vector<unsigned> validRefAlignment;
    std::vector<unsigned> validSnpAlignment;

    // Process reference and variant genome mapping
    std::cout << "Process reference off-targets" << std::endl;
    readBamFile(refOffTargets, refBamPath, faiIndexRef, numMismatches);
    sortSnpRegionsByChr(sortedIndexAllChr, chrMap, snpInfoTable, threads);
    filterRefAlignment(validRefAlignment, sortedIndexAllChr, chrMap, snpInfoTable, refOffTargets, onTargets, seqLength, threads);

    std::cout << "Process variant off-targets" << std::endl;
    readBamFile(snpOffTargets, snpBamPath, faiIndexSnp, numMismatches);
    filterSnpAlignment(validSnpAlignment, snpOffTargets, onTargets, snpInfoTable, seqLength);

    // Read TUSCAN activity file in
    std::map<CharString, double> onTargetActivity;
    readTuscanResult(onTargetActivity, tuscanPath);

    // Write combined results to a file
    std::ofstream mergedResults;
    mergedResults.open(outputPath);

    std::ofstream featureMatrix;
    featureMatrix.open(featureMatrixPath);

    if (mergedResults.is_open())
    {
        // Write column names output
        mergedResults << "#Chr\t" << "Start\t" << "End\t" << "Targetsite\t" << "Score\t" << "Strand\t" << "Sequence\t" << "Mismatch_Number\t" << "Mismatch_Positions\t" << "Variants\n";

        // Write column names feature matrix
        std::vector<std::string> featureNames;
        getFeatureNames(featureNames, seqLength);

        for (unsigned i = 0; i < (featureNames.size() - 1); i++)
        {
            featureMatrix << featureNames[i] << "\t";
        }
        featureMatrix << featureNames[featureNames.size() - 1] << "\n";

        // Write reference results
        std::vector<int> perfectMatch{-1};
        double currOntargetActivity;
        std::string offTargetName;

        for (unsigned i = 0; i < validRefAlignment.size(); i++)
        {
            // Count off-targets and get current off-target count for name
            offTargetCount.at(refOffTargets[validRefAlignment[i]].target)++;
            offTargetName = (std::string) toCString(refOffTargets[validRefAlignment[i]].target) + '_' + std::to_string(offTargetCount.at(refOffTargets[validRefAlignment[i]].target));

            mergedResults << refOffTargets[validRefAlignment[i]].chr << "\t";
            mergedResults << refOffTargets[validRefAlignment[i]].pos << "\t";
            mergedResults << (refOffTargets[validRefAlignment[i]].pos + 23) << "\t";
            mergedResults << offTargetName << "\t";
            mergedResults << ".\t";
            mergedResults << refOffTargets[validRefAlignment[i]].strand << "\t";
            mergedResults << refOffTargets[validRefAlignment[i]].sequence << "\t";

            if (refOffTargets[validRefAlignment[i]].mismatchPos == perfectMatch)
            {
                mergedResults << 0 << "\t\t";
            }
            else
            {
                mergedResults << refOffTargets[validRefAlignment[i]].mismatchPos.size() << "\t";
                for (unsigned j = 0; j < (refOffTargets[validRefAlignment[i]].mismatchPos.size() - 1); j++)
                {
                    mergedResults << refOffTargets[validRefAlignment[i]].mismatchPos[j] << ',';
                }
                mergedResults << refOffTargets[validRefAlignment[i]].mismatchPos[refOffTargets[validRefAlignment[i]].mismatchPos.size() - 1] << "\t";
            }
            mergedResults << refOffTargets[validRefAlignment[i]].snpType << "\n";

            // Feature matrix record
            std::vector<unsigned> features;
            featureMatrixRecord(features, onTargets[refOffTargets[validRefAlignment[i]].target].sequence, refOffTargets[validRefAlignment[i]].sequence);

            // Row name
            featureMatrix << offTargetName << "\t";

            // All sequence features
            for (unsigned j = 0; j < features.size(); j++)
            {
                featureMatrix << features[j] << "\t";
            }

            // On-target activity
            featureMatrix << onTargetActivity.at(refOffTargets[validRefAlignment[i]].target) << "\n";
        }

        // Write variant results
        for (unsigned i = 0; i < validSnpAlignment.size(); i++)
        {
            // Count off-targets and get current off-target count for name
            offTargetCount.at(snpOffTargets[validSnpAlignment[i]].target)++;
            offTargetName = (std::string) toCString(snpOffTargets[validSnpAlignment[i]].target) + '_' + std::to_string(offTargetCount.at(snpOffTargets[validSnpAlignment[i]].target));

            mergedResults << snpOffTargets[validSnpAlignment[i]].chr << "\t";
            mergedResults << snpOffTargets[validSnpAlignment[i]].pos << "\t";
            mergedResults << (snpOffTargets[validSnpAlignment[i]].pos + 23) << "\t";
            mergedResults << offTargetName << "\t";
            mergedResults << ".\t";
            mergedResults << snpOffTargets[validSnpAlignment[i]].strand << "\t";
            mergedResults << snpOffTargets[validSnpAlignment[i]].sequence << "\t";

            if (snpOffTargets[validSnpAlignment[i]].mismatchPos == perfectMatch)
            {
                mergedResults << 0 << "\t\t";
            }
            else
            {
                mergedResults << snpOffTargets[validSnpAlignment[i]].mismatchPos.size() << "\t";
                for (unsigned j = 0; j < (snpOffTargets[validSnpAlignment[i]].mismatchPos.size() - 1); j++)
                {
                    mergedResults << snpOffTargets[validSnpAlignment[i]].mismatchPos[j] << ',';
                }
                mergedResults << snpOffTargets[validSnpAlignment[i]].mismatchPos[snpOffTargets[validSnpAlignment[i]].mismatchPos.size() - 1] << "\t";
            }

            mergedResults << snpOffTargets[validSnpAlignment[i]].snpType << "\n";

            // Feature matrix record
            std::vector<unsigned> features;
            featureMatrixRecord(features, onTargets[snpOffTargets[validSnpAlignment[i]].target].sequence, snpOffTargets[validSnpAlignment[i]].sequence);

            // Row name
            featureMatrix << offTargetName << "\t";

            // All sequence features
            for (unsigned j = 0; j < features.size(); j++)
            {
                featureMatrix << features[j] << "\t";
            }

            // On-target activity
            featureMatrix << onTargetActivity.at(snpOffTargets[validSnpAlignment[i]].target) << "\n";
        }

        mergedResults.close();
        featureMatrix.close();
        std::cout << "Merging output files finished" << std::endl;
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open output file.");
    }
}


/*!
 * @fn processRefOnly
 * @brief Processes reference only mapping and computes MIT score for output
 *
 * @signature void mergeResults(outputPath, featureMatrixPath, refBamPath, ontargetPath, indexPathRef, tuscanPath,
 *                              numMismatches, seqLength)
 *
 * @param[in]       outputPath          Path of the merged output file
 * @param[in]       featureMatrixPath   Path of the merged output file
 * @param[in]       refBamPath          Path to reference genome BAM/SAM file
 * @param[in]       ontargetPath        Path to ontarget fasta file
 * @param[in]       indexPathRef        Path for the genome, index should be built of or index exists in same directory
 * @param[in]       tuscanPath          Path to on-target activity (TUSCAN regression) file
 * @param[in]       numMismatches       Maximum number of mismatches allowed during preceeding read mapping
 * @param[in]       seqLength           Length of the target sequences that were searched
 *
 * This function merges the two output files from mapping against reference and SNP genome. In addition, it is secured
 * that no duplicates, matches to N region in combined SNP genome, matches to the reference in SNP regions or
 * off-targets with PAM mismatches are included into the final output file.
 * The final output contains for each potential off-target the bed6 records where name is the on-target name, the
 * sequence, the positions of mismatches (if existing) and specific tags if variants are included in the sequence.
 */
void processRefOnly(char const * outputPath, char const * refBamPath, char const * ontargetPath,
char const * indexPathRef, char const * tuscanPath, unsigned const & numMismatches, unsigned const & seqLength)
{
    // Open fai index reference
    FaiIndex faiIndexRef;
    if (!open(faiIndexRef, indexPathRef))
    {
        if (!build(faiIndexRef, indexPathRef))
        {
            throw std::runtime_error("ERROR: Reference index could not be loaded or built.");
        }
        if (!save(faiIndexRef))
        {
            throw std::runtime_error("ERROR: Variant index could not be written do disk.");
        }
    }

    // Read in all potential off-targets from reference and SNP genome
    std::vector<PotentialOffTarget> refOffTargets;

    std::cout << "Read reference BAM file" << std::endl;
    readBamFile(refOffTargets, refBamPath, faiIndexRef, numMismatches);

    // Read in on-targets and prepare counting number of off-targets per on-target
    std::map<CharString, PotentialOffTarget> onTargets;
    std::map<CharString, unsigned> offTargetCount;

    readOntargets(onTargets, offTargetCount, ontargetPath, faiIndexRef, numMismatches);

    // Read TUSCAN activity file in
    std::map<CharString, double> onTargetActivity;
    readTuscanResult(onTargetActivity, tuscanPath);

    // Write results to a file
    std::ofstream refResults;
    refResults.open(outputPath);

    bool valid;
    if (refResults.is_open())
    {
        refResults << "#Chr\t" << "Start\t" << "End\t" << "Targetsite\t" << "Score\t" << "Strand\t" << "Sequence\t" << "Mismatch_Number\t" << "Mismatch_Positions\n";

        std::vector<int> perfectMatch{-1};
        double currOntargetActivity;
        std::string offTargetName;

        for (unsigned i = 0; i < refOffTargets.size(); i++)
        {
            valid = true;
            if (comp(refOffTargets[i], onTargets.at(refOffTargets[i].target)))
            {
                valid = false;
            }

            if (valid)
            {
                // Count off-targets and get current off-target count for name
                offTargetCount.at(refOffTargets[i].target)++;
                offTargetName = (std::string) toCString(refOffTargets[i].target) + '_' + std::to_string(offTargetCount.at(refOffTargets[i].target));

                refResults << refOffTargets[i].chr << "\t";
                refResults << refOffTargets[i].pos << "\t";
                refResults << (refOffTargets[i].pos + 23) << "\t";
                refResults << offTargetName << "\t";
                refResults << calcMitScore(refOffTargets[i].mismatchPos) << "\t";
                refResults << refOffTargets[i].strand << "\t";
                refResults << refOffTargets[i].sequence << "\t";

                if (refOffTargets[i].mismatchPos == perfectMatch)
                {
                    refResults << 0 << "\t\n";
                }
                else
                {
                    refResults << refOffTargets[i].mismatchPos.size() << "\t";
                    for (unsigned j = 0; j < (refOffTargets[i].mismatchPos.size() - 1); j++)
                    {
                        refResults << refOffTargets[i].mismatchPos[j] << ',';
                    }
                    refResults << refOffTargets[i].mismatchPos[refOffTargets[i].mismatchPos.size() - 1] << "\n";
                }
            }
        }

        refResults.close();
        std::cout << "Writing reference output finished." << std::endl;
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open output file.");
    }
}

/*!
 * @fn processRefOnly
 * @brief Processes reference only mapping and writes feature matrix as second output
 *
 * @signature void mergeResults(outputPath, featureMatrixPath, refBamPath, ontargetPath, indexPathRef, tuscanPath,
 *                              numMismatches, seqLength)
 *
 * @param[in]       outputPath          Path of the merged output file
 * @param[in]       featureMatrixPath   Path of the merged output file
 * @param[in]       refBamPath          Path to reference genome BAM/SAM file
 * @param[in]       ontargetPath        Path to ontarget fasta file
 * @param[in]       indexPathRef        Path for the genome, index should be built of or index exists in same directory
 * @param[in]       tuscanPath          Path to on-target activity (TUSCAN regression) file
 * @param[in]       numMismatches       Maximum number of mismatches allowed during preceeding read mapping
 * @param[in]       seqLength           Length of the target sequences that were searched
 *
 * This function merges the two output files from mapping against reference and SNP genome. In addition, it is secured
 * that no duplicates, matches to N region in combined SNP genome, matches to the reference in SNP regions or
 * off-targets with PAM mismatches are included into the final output file.
 * The final output contains for each potential off-target the bed6 records where name is the on-target name, the
 * sequence, the positions of mismatches (if existing) and specific tags if variants are included in the sequence.
 */
void processRefOnly(char const * outputPath, char const * featureMatrixPath, char const * refBamPath,
char const * ontargetPath, char const * indexPathRef, char const * tuscanPath, unsigned const & numMismatches,
unsigned const & seqLength)
{
    // Open fai index reference
    FaiIndex faiIndexRef;
    if (!open(faiIndexRef, indexPathRef))
    {
        if (!build(faiIndexRef, indexPathRef))
        {
            throw std::runtime_error("ERROR: Reference index could not be loaded or built.");
        }
        if (!save(faiIndexRef))
        {
            throw std::runtime_error("ERROR: Variant index could not be written do disk.");
        }
    }

    // Read in all potential off-targets from reference and SNP genome
    std::vector<PotentialOffTarget> refOffTargets;

    std::cout << "Read reference BAM file" << std::endl;
    readBamFile(refOffTargets, refBamPath, faiIndexRef, numMismatches);

    // Read in on-targets and prepare counting number of off-targets per on-target
    std::map<CharString, PotentialOffTarget> onTargets;
    std::map<CharString, unsigned> offTargetCount;

    readOntargets(onTargets, offTargetCount, ontargetPath, faiIndexRef, numMismatches);

    // Read TUSCAN activity file in
    std::map<CharString, double> onTargetActivity;
    readTuscanResult(onTargetActivity, tuscanPath);

    // Write results to a file
    std::ofstream refResults;
    refResults.open(outputPath);

    std::ofstream featureMatrix;
    featureMatrix.open(featureMatrixPath);

    bool valid;
    if (refResults.is_open() && featureMatrix.is_open())
    {
        refResults << "#Chr\t" << "Start\t" << "End\t" << "Targetsite\t" << "Score\t" << "Strand\t" << "Sequence\t" << "Mismatch_Number\t" << "Mismatch_Positions\n";
        std::vector<std::string> featureNames;
        getFeatureNames(featureNames, seqLength);

        for (unsigned i = 0; i < (featureNames.size() - 1); i++)
        {
            featureMatrix << featureNames[i] << "\t";
        }
        featureMatrix << featureNames[featureNames.size() - 1] << "\n";

        std::vector<int> perfectMatch{-1};
        double currOntargetActivity;
        std::string offTargetName;

        for (unsigned i = 0; i < refOffTargets.size(); i++)
        {
            valid = true;
            if (comp(refOffTargets[i], onTargets.at(refOffTargets[i].target)))
            {
                valid = false;
            }

            if (valid)
            {
                // Count off-targets and get current off-target count for name
                offTargetCount.at(refOffTargets[i].target)++;
                offTargetName = (std::string) toCString(refOffTargets[i].target) + '_' + std::to_string(offTargetCount.at(refOffTargets[i].target));

                refResults << refOffTargets[i].chr << "\t";
                refResults << refOffTargets[i].pos << "\t";
                refResults << (refOffTargets[i].pos + 23) << "\t";
                refResults << offTargetName << "\t";
                refResults << ".\t";
                refResults << refOffTargets[i].strand << "\t";
                refResults << refOffTargets[i].sequence << "\t";

                if (refOffTargets[i].mismatchPos == perfectMatch)
                {
                    refResults << 0 << "\t\n";
                }
                else
                {
                    refResults << refOffTargets[i].mismatchPos.size() << "\t";
                    for (unsigned j = 0; j < (refOffTargets[i].mismatchPos.size() - 1); j++)
                    {
                        refResults << refOffTargets[i].mismatchPos[j] << ',';
                    }
                    refResults << refOffTargets[i].mismatchPos[refOffTargets[i].mismatchPos.size() - 1] << "\n";
                }

                // Feature matrix record
                std::vector<unsigned> features;
                featureMatrixRecord(features, onTargets[refOffTargets[i].target].sequence, refOffTargets[i].sequence);

                // Row name
                featureMatrix << offTargetName << "\t";

                // All sequence features
                for (unsigned j = 0; j < features.size(); j++)
                {
                    featureMatrix << features[j] << "\t";
                }

                // On-target activity
                featureMatrix << onTargetActivity.at(refOffTargets[i].target) << "\n";
            }
        }
        refResults.close();
        featureMatrix.close();

        std::cout << "Writing reference output finished." << std::endl;
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open output file.");
    }
}

}
