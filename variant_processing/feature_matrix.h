// ============================================================================
// Compute feature matrix record
// ============================================================================

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace seqan
{

/*!
 * @fn featureMatrixRecord
 * @brief Computes feature matrix record for one off-target and the corresponding on-target
 *
 * @signature void featureMatrixRecord(features, onTarget, offTarget)
 *
 * @param[in,out]   features            Vector containing all features
 * @param[in]       onTarget            On-target
 * @param[in]       offTarget           Off-target
 */
void featureMatrixRecord(std::vector<unsigned> & features, Dna5String const & onTarget, Dna5String const & offTarget)
{
    // On-target activity is added later
    // Position 0         = NM                             (1 feature)
    // Position 1   - 21  = Mismatch positions             (21 features)
    // Position 22  - 33  = Mismatch type                  (12 features)
    // Position 34        = Transition number              (1 feature)
    // Position 35        = Transversion number            (1 feature)
    // Position 36  - 115 = Single DNA letters             (80 features)
    // Position 116 - 119 = PAM letters                    (4 features)
    // Position 120 - 423 = Paired DNA letters             (304 features)
    // Position 424 - 439 = Number of paired DNA letters   (16 features)
    // Position 440       = Adjacent mismatch number       (1 feature)
    // Position 441       = Seed mismatch number           (1 feature)

    features.resize(442, 0);

    std::map<Dna5String, unsigned> dnaPairs = {{"AA", 0}, {"AC", 1}, {"AG", 2}, {"AT", 3}, {"CA", 4}, {"CC", 5}, {"CG", 6},
                                               {"CT", 7}, {"GA", 8}, {"GC", 9}, {"GG", 10}, {"GT", 11}, {"TA", 12}, {"TC", 13},
                                               {"TG", 14}, {"TT", 15}};
    std::map<Dna5String, unsigned> mismatchTypes = {{"AC", 0}, {"AG", 1}, {"AT", 2}, {"CA", 3}, {"CG", 4}, {"CT", 5},
                                                    {"GA", 6}, {"GC", 7}, {"GT", 8}, {"TA", 9}, {"TC", 10}, {"TG", 11}};
    std::set<Dna5String> transitions = {"AG", "CT", "GA", "TC"};

    // Iterator to loop over string once
    bool precMismatch = false;
    Dna5String currMismatch, currPair;

    for (unsigned i = 0; i < (length(offTarget) - 2); ++i)
    {
        // Paired letters
        if (i < 19)
        {
            currPair = infix(offTarget, i, i + 2);
            features[120 + i * dnaPairs.size() + dnaPairs[currPair]] = 1;
            features[424 + dnaPairs[currPair]]++;
        }

        // Single letters
        switch ((char) offTarget[i])
        {
            case 'A':
                // 4 is alphabet size
                // Dna5String needs to be chosen anyway as Ns might exist and filtered out earlier but types should be consistent later
                features[36 + i * 4] = 1;
                break;
            case 'C':
                features[36 + i * 4 + 1] = 1;
                break;
            case 'G':
                features[36 + i * 4 + 2] = 1;
                break;
            case 'T':
                features[36 + i * 4 + 3] = 1;
                break;
            default:
                features[36 + i * 4] = 1;
                break;
        }

        // Mismatches
        if (onTarget[i] != offTarget[i])
        {
            // Total number of mismatches
            features[0]++;

            // Mismatch positions
            features[i + 1] = 1;

            if (i > 7 && i < 20)
            {
                // Seed mismatches
                features[441]++;
            }

            if (precMismatch)
            {
                features[440]++;
            }

            precMismatch = true;

            currMismatch = onTarget[i];
            append(currMismatch, offTarget[i]);

            // Transition and transversion number
            if (transitions.find(currMismatch) != transitions.end())
            {
                features[34]++;
            }
            else
            {
                features[35]++;
            }
            features[22 + mismatchTypes[currMismatch]] = 1;
        }
        else
        {
            precMismatch = false;
        }
    }
}

/*!
 * @fn getFeatureNames
 * @brief Feature names for feature matrix columns.
 *
 * @signature void getFeatureNames()
 *
 * @param[in,out]           featureNames        Vector containing the feature names
 * @param[in]               seqLength           Length of the off-targets
 *
 * This needs to be adapted if ever new features should be added.
 * In addition, the featureMatrixRecord function needs to be changed as exact positions are predefined.
 */
void getFeatureNames(std::vector<std::string> & featureNames, unsigned seqLength)
{
    // Position 0         = NM                             (1 feature)
    // Position 1   - 21  = Mismatch positions             (21 features)
    // Position 22  - 33  = Mismatch type                  (12 features)
    // Position 34        = Transition number              (1 feature)
    // Position 35        = Transversion number            (1 feature)
    // Position 36  - 115 = Single DNA letters             (80 features)
    // Position 116 - 119 = PAM letters                    (4 features)
    // Position 120 - 423 = Paired DNA letters             (304 features)
    // Position 425 - 439 = Number of paired DNA letters   (16 features)
    // Position 440       = Adjacent mismatch number       (1 feature)
    // Position 441       = Seed mismatch number           (1 feature)
    // Position 442       = On-target activity             (1 feature)

    featureNames.resize(443);
    std::vector<std::string> const mismatchTypes = {"AtoC", "AtoG", "AtoT", "CtoA", "CtoG", "CtoT", "GtoA", "GtoC", "GtoT", "TtoA", "TtoC", "TtoG"};
    std::vector<std::string> const dnaLetters = {"A", "C", "G", "T"};
    std::vector<std::string> const dnaPairs = {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"};

    featureNames[0] = "totalMismatches";

    for (unsigned i = 1; i < (seqLength - 1); i++)
    {
        featureNames[i] = "mismatchPos" + std::to_string(i);
    }

    for (unsigned i = 22; i < (22 + mismatchTypes.size()); i++)
    {
        featureNames[i] = mismatchTypes[i - 22];
    }

    featureNames[34] = "transitionNumber";
    featureNames[35] = "transversionNumber";

    for (unsigned i = 1; i < (seqLength - 2); i++)
    {
        for (unsigned j = 0; j < dnaLetters.size(); j++)
        {
            featureNames[36 + (i - 1) * dnaLetters.size() + j] = dnaLetters[j] + std::to_string(i);
        }
    }

    featureNames[116] = "PAMA";
    featureNames[117] = "PAMC";
    featureNames[118] = "PAMG";
    featureNames[119] = "PAMT";

    for (unsigned i = 1; i < (seqLength - 3); i++)
    {
        for (unsigned j = 0; j < dnaPairs.size(); j++)
        {
            featureNames[120 + (i - 1) * dnaPairs.size() + j] = dnaPairs[j] + std::to_string(i);
        }
    }

    for (unsigned i = 424; i < (424 + dnaPairs.size()); i++)
    {
        featureNames[i] = dnaPairs[i - 424];
    }

    featureNames[440] = "adjacentMismatches";
    featureNames[441] = "seedMismatches";
    featureNames[442] = "ontargetActivity";
}

void readTuscanResult(std::map<CharString, double> & onTargetActivity, char const * inputPath)
{
    std::string record;
    std::ifstream in(inputPath);

    std::string target, sequence;
    double score;

    if (in.is_open())
    {
        while (getline(in, record))
        {
            std::istringstream is(record);
            if (is >> target >> sequence >> score)
            {
                onTargetActivity.insert(std::make_pair((CharString) target, score));
            }
        }
        in.close();
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open on-target activity file.");
    }
}

}
