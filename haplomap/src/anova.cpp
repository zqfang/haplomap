//
// Created by Zhuoqing Fang on 7/1/20.
//
#include "ghmap.h"
#include "gsl/gsl_cdf.h"


std::vector<float> sumVector(std::vector<float> &vec) {
    std::vector<float> res(1,0.0F);
    for (float & v: vec )
        res[0] += v;
    return res;
}

void subVector(std::vector<float> &vec, float value) {
    //std::vector<float> res;
    for (float & v: vec ){
        v -= value;
        //res.push_back(v);
    }
    //return res;
}

void ANOVA(vector<vector<float>> &phenvec, char *pattern, float &FStat, float &pvalue, float &effect)
{
    int numHaplo = numHaplotypes(pattern);
    // array haplotype -> num strains in haplotype.
    vector<int> haploNum(numHaplo, 0);

    // array haplotype -> mean (vector<float>) for each haplotype
    vector<vector<float>> haploMean(numHaplo, vector<float>(numCategories, 0.0F)); // size 1

    float numDefined = 0.0F; // numbers of strains without ?
    float numDefinedDataPoints = 0.0F; // number of all individual data point without '?' strain
    vector<float> sumDefined(numCategories, 0.0F); //size 1

    std::vector<std::vector<float>> sumStrains;

    // Compute haplotype means
    for (int str1 = 0; str1 < numStrains; str1++)
    {
        int hap = pattern[str1]; // 0,1,2,3,4
        vector<float> &phen = phenvec[str1]; // multivalue?
        if ('?' != hap)
        {
            numDefined ++;
            numDefinedDataPoints += phen.size();
            sumStrains.push_back(sumVector(phen));
            //addVectors(sumDefined, sumStrains.back());
            //haploNum[hap]++;

            //addVectors(haploMean[hap], phen); // temporarily, the total, not mean.
            addVectors(haploMean[hap], sumStrains.back());
            haploNum[hap] += phen.size();
        }
    }

    for (int hap = 0; hap < numHaplo; hap++)
    {
        scaleVector(haploMean[hap], 1.0 / haploNum[hap]); // get the mean
    }

    if (traceFStat)
    {
        cout << "numDefined = " << numDefined << ", sumDefined = " << sumDefined.back() << endl;
        cout << "haploNum[] = [";
        for (int hap = 0; hap < numHaplo; hap++)
        {
            cout << haploNum[hap] << " ";
        }
        cout << "]" << endl;

        cout << "haploMean[] = [";
        for (int hap = 0; hap < numHaplo; hap++)
        {
            cout << haploMean[hap].back() << " ";
        }
        cout << "]" << endl;
    }

    float SSW = 0.0F;
    for (int str1 = 0; str1 < numStrains; str1++)
    {
        int hap = pattern[str1];
        if ('?' != hap)
        {
            vector<float> resid = phenvec[str1];
            subVector(resid, haploMean[hap].back());
            //subtractVectors(resid, haploMean[hap]);
            float tmpdot = dotVectors(resid, resid);
            SSW += tmpdot;
        }
    }

    // SSB -- between sum of squares (sum over haplotypes hapsize*(hapmean-mean)^2
    float SSB = 0.0;
    vector<float> mean(1, 0.0F); // mean of all data
    for (auto & s: sumStrains)
        mean[0] += s.back();
    // scaleVector(mean, 1.0 / numDefined);
    scaleVector(mean, 1.0 / numDefinedDataPoints);

    for (int hap = 0; hap < numHaplo; hap++)
    {
        vector<float> diff = haploMean[hap]; // copy so we don't destroy haploMeans
        subtractVectors(diff, mean);         // (haplotype mean) - mean
        // FIXME: this is a wrong formula? why haploNum[hap] here?
        float sq = haploNum[hap] * dotVectors(diff, diff);
        // Should just sq = dotVectors(diff, diff);
        //float sq = dotVectors(diff, diff);
        SSB += sq;
    }
    // FIXME:
    // my quick and dirty hack to penalize missing alleles.
    // simulates additional error for each missing value (but ignores
    // degrees of freedom).
    SSW += 1.1 * (numStrains - numDefined);

    // mean square within
    // df within is (numDefined-numHaplo)
    float dfW = numDefinedDataPoints - numHaplo; // degrees of freedom within
    float dfB = numHaplo - 1;          // degrees of freedom between.
    if (dfW == 0.0)
    {
        // This happens when numDefined = numHaplo, which occurs rarely when there are
        // lots of undefined strains and lots of haplotypes.
        // This causes errors in gsl, and I don't know what the right thing to do is,
        // so just punt.  We won't get a match for this.
        pvalue = 1.0;
        effect = 0.0;
        return;
    }
    float MSW = SSW / dfW;
    float MSB = SSB / dfB;

    // This formula is the same as Peltz.  So, I think I've seen two totally
    // different formulas for omega^2
    // WARNING: This divides by 0 if SSW is 0.
    // Which seems to work ok (F <- "inf").
    FStat = MSB / MSW;

    // out parameter for pvalue
    pvalue = (float)gsl_cdf_fdist_Q((double)FStat,
                                    (double)(numHaplo - 1),
                                    (double)(numDefined - numHaplo));

    // Genetic effect
    // Oh wow!  omega^2 is the "coefficient of determination"!
    // http://faculty.chass.ncsu.edu/garson/PA765/anova.htm#anova2
    effect = (float)((SSB - (numHaplo - 1) * MSW) / (SSW + SSB + MSW));

    if (traceFStat)
    {
        cout << "SSW = " << SSW
             << ", SSB = " << SSB
             << ", MSW = " << MSW
             << ", MSB = " << MSB
             << ", FStat = " << FStat
             << ", pval = " << pvalue
             << ", genetic effect = " << effect
             << endl;
    }
}