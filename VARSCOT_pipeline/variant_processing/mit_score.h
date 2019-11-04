// ============================================================================
// MIT score implementation from 'scores of single hits' on http://crispr.mit.edu/about
// ============================================================================

#pragma once

#include <cmath>
#include <vector>

#include <seqan/sequence.h>

double calcMitScore(std::vector<int> const & mismatchPos)
{
    double s, s1, s2, s3, avgDist;

    std::vector<int> const perfectMatch{-1};
    unsigned nm;

    if (mismatchPos == perfectMatch)
    {
        nm = 0;
    }
    else
    {
        // Exclude mismatch in the PAM
        if (mismatchPos[mismatchPos.size() - 1] < 20)
        {
            nm = mismatchPos.size();
        }
        else
        {
            nm = mismatchPos.size() - 1;
        }

        s3 = (double) 1 / (double) std::pow(nm, 2);
    }

    if(nm==0)
    {
        return 100;
    }
    std::vector<double> const matrixM{0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583};

    s1 = 1;
    std::vector<unsigned> dist;
    dist.reserve(nm - 1);

    for (unsigned i = 0; i < nm; i++)
    {
        s1 *= (1 - matrixM[mismatchPos[i]]);
        if (i > 0)
        {
            dist.push_back(mismatchPos[i] - mismatchPos[i - 1]);
        }
    }

    if (nm < 2)
    {
        s2 = 1;
    }
    else
    {
        avgDist = (double) std::accumulate(dist.begin(), dist.end(), 0) / (double) dist.size();
        s2 = 1 / (((19 - avgDist) / 19) * 4 + 1);
    }
    s = s1 * s2 * s3 * 100;
    return s;
}

