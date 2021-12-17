//
// Created by Zhuoqing Fang on 7/1/20.
//

// #include <numeric>
// #include "gsl/gsl_cdf.h"
// #include "gsl/gsl_math.h"
// #include "ghmap.h"

//// old code
//void ANOVA(std::vector<std::vector<float>> &phenvec, char *pattern, float &FStat, float &pvalue, float &effect)
//{
//    int numHaplo = numHaplotypes(pattern);
//    // array haplotype -> num strains in haplotype.
//    vector<int> haploNum(numHaplo, 0);
//
//    // array haplotype -> mean (vector<float>) for each haplotype
//    vector<vector<float>> haploMean(numHaplo, vector<float>(numCategories, 0.0F)); // size 1
//
//    float numDefined = 0.0F; // number of all individual data point
//    vector<float> sumDefined(numCategories, 0.0F); //size 1
//
//    // Compute haplotype means
//    for (int str1 = 0; str1 < numStrains; str1++)
//    {
//        int hap = pattern[str1]; // 0,1,2,3,4
//        vector<float> &phen = phenvec[str1]; // multivalue?
//        if ('?' != hap)
//        {
//            numDefined += phen.size();
//            addVectors(sumDefined, phen);
//            haploNum[hap]++;
//
//            addVectors(haploMean[hap], phen); // temporarily, the total, not mean.
//        }
//    }
//
//    for (int hap = 0; hap < numHaplo; hap++)
//    {
//        scaleVector(haploMean[hap], 1.0 / haploNum[hap]); // get the mean
//    }
//
//    if (traceFStat)
//    {
//        // FIXME: sumDefined
//        cout << "numDefined = " << numDefined << ", sumDefined = " << /*sumDefined <<*/ endl;
//        cout << "haploNum[] = [";
//        for (int hap = 0; hap < numHaplo; hap++)
//        {
//            cout << haploNum[hap] << " ";
//        }
//        cout << "]" << endl;
//
//        cout << "haploMean[] = [";
//        for (int hap = 0; hap < numHaplo; hap++)
//        {// FIXME:
//            //cout << haploMean[hap] << " ";
//        }
//        cout << "]" << endl;
//    }
//
//    float SSW = 0.0F;
//    for (int str1 = 0; str1 < numStrains; str1++)
//    {
//        int hap = pattern[str1];
//        if ('?' != hap)
//        {
//            vector<float> resid = phenvec[str1];
//            subtractVectors(resid, haploMean[hap]);
//            float tmpdot = dotVectors(resid, resid);
//
//            SSW += tmpdot;
//        }
//    }
//
//    // SSB -- between sum of squares (sum over haplotypes hapsize*(hapmean-mean)^2
//    float SSB = 0.0;
//    vector<float> &mean = sumDefined; // sumDefined will be the mean of all values.
//    scaleVector(mean, 1.0 / numDefined);
//
//    for (int hap = 0; hap < numHaplo; hap++)
//    {
//        vector<float> diff = haploMean[hap]; // copy so we don't destroy haploMeans
//        subtractVectors(diff, mean);         // (haplotype mean) - mean
//        float sq = haploNum[hap] * dotVectors(diff, diff);
//        SSB += sq;
//    }
//
//    // my quick and dirty hack to penalize missing alleles.
//    // simulates additional error for each missing value (but ignores
//    // degrees of freedom).
//    SSW += 1.1 * (numStrains - numDefined);
//
//    // mean square within
//    // df within is (numDefined-numHaplo)
//    float dfW = numDefined - numHaplo; // degrees of freedom within
//    float dfB = numHaplo - 1;          // degrees of freedom between.
//    if (dfW == 0.0)
//    {
//        // This happens when numDefined = numHaplo, which occurs rarely when there are
//        // lots of undefined strains and lots of haplotypes.
//        // This causes errors in gsl, and I don't know what the right thing to do is,
//        // so just punt.  We won't get a match for this.
//        pvalue = 1.0;
//        effect = 0.0;
//        return;
//    }
//    float MSW = SSW / dfW;
//    float MSB = SSB / dfB;
//
//    // This formula is the same as Peltz.  So, I think I've seen two totally
//    // different formulas for omega^2
//    // WARNING: This divides by 0 if SSW is 0.
//    // Which seems to work ok (F <- "inf").
//    FStat = MSB / MSW;
//    // out parameter for pvalue
//    pvalue = (float)gsl_cdf_fdist_Q((double)FStat,
//                                    (double)(numHaplo - 1),
//                                    (double)(numDefined - numHaplo));
//
//    // Genetic effect
//    // Oh wow!  omega^2 is the "coefficient of determination"!
//    // http://faculty.chass.ncsu.edu/garson/PA765/anova.htm#anova2
//    effect = (float)((SSB - (numHaplo - 1) * MSW) / (SSW + SSB + MSW));
//
//    if (traceFStat)
//    {
//        cout << "SSW = " << SSW
//             << ", SSB = " << SSB
//             << ", MSW = " << MSW
//             << ", MSB = " << MSB
//             << ", FStat = " << FStat
//             << ", pval = " << pvalue
//             << ", genetic effect = " << effect
//             << endl;
//    }
// }
// SST = SSW + SSB
// (SSB/(k-1))/(SSW/(n-k)) -- k = numhaplotypes, n == numstrains (defined?)

// Trying to correct between variations in notation/typos:  AP stats notes say:
// F0 = MST(Between)/MSE(within)
// MST(Between) = SST(between)/(df(Between)
// MST(Within) = SST(Within)/(df(within))
// Ok, so MSE(Between) = SSB/(k-1)
// Ok, so MSE(Within) = SSW/(n-k)

// "effect" "treatment" "between" seem to be roughly the same.
// df(effect) is k-1
// "error" "within" are equivalent.
// df(error) is n-k

// MSE - mean square error SSE/df(error)

// omega^2 = (SSE k- df(effect)(MSE)) / (MSE + SST)
// Peltz: omega^2 = (SSB - (k-1)*MSE)/(SST + MSE)
//  INCONSISTENCY: I thought SSE was SSW