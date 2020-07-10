//
// Created by Zhuoqing Fang on 6/16/20.
//

#include <string>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "ghmap.h"
#include "manova.h"
#include "fdr.h"
#include "eigen.h"

using namespace std;

template <typename NumericType>
void print_vec(std::vector<NumericType> &vec) {
    int i = 0;
    for (auto & v: vec){
        std::cout<<v<<" ";
        i++;
        if (i % 10 == 0)
            std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

// debugging
void print_gslvec(gsl_vector* M)
{
    std::cout<<"trace vector: "<<std::endl;
    for (int i = 0; i < M->size; ++i)
    {
        std::cout<<gsl_vector_get(M, i)<<" ";
    }
    std::cout<<std::endl;
}

void print_gslmat(gsl_matrix* M){
    std::cout << "trace matrix: "<<std::endl;
    for (int i = 0; i < M->size1; ++i){
        for (int j=0; j < M->size2; ++j) {
            std::cout<< gsl_matrix_get(M, i, j) << " ";
        }
        std::cout<<std::endl;
    }
}




TEST(READ_PHENO, DISABLED_cat_and_quant_read) {
    int numStrains = 12;
    vector<vector<float>> phenvec(100);
    bool isCategorical = false;
    char* path = (char *)"../../data/test.trait.txt";
    if (isCategorical)
    {
        readCPhenotypes(path, phenvec);
    }
    else
    {
        readQPhenotypes(path, phenvec);
    }

    cout << "Phenotype vectors for strains" << endl;
    for (int str = 0; str < phenvec.size(); str++)
    {
        for (int c=0; c < phenvec[0].size(); c++)
            cout << phenvec[str][c] << " ";
        cout<<endl;
    }
    cout << endl;

}

TEST(READ_BLOCK, read_block_summary) {
    char* path = (char *)"../../data/MPD_39504/chrX.hblocks.txt";
    std::vector<BlockSummary *> xblocks;
    readBlockSummary(path, NULL, false);
    BlockSummary *pBlock = xblocks.back();

    EXPECT_NE(pBlock, nullptr);

    cout << "Block: " << pBlock->chrName << "\t" << pBlock->blockIdx << "\t"<<endl;
    showPattern(pBlock->pattern);
}




TEST(FDR_RUN, DISABLED_fdr_test){
    vector<double> pval = { 0.001063732,0.003773763,0.006233452,0.006535265,0.01132544,
                          0.01838011,0.01936604,0.02866878,0.03798249,0.03986961,
                          0.04266445,0.04405806,0.05517378,0.05996954,0.07529589,
                          0.08391378,0.08796254,0.09734611,0.1029043,
                          0.1111793,0.1115418,0.1193652,0.1869531,0.1902924,0.3105688 };
    std::vector<double> padj(pval.size(), 1.0);
    std::vector<bool> reject(pval.size(), false);
    reject = bh_fdr(pval, padj, 0.05);
    std::cout<<"Pvalue:"<<std::endl;
    print_vec(pval);
    std::cout<<"FDR"<<std::endl;
    print_vec(padj);
}

TEST(CorMat_READ, DISABLED_cormat_test)
{
    //EXPECT_TRUE(aov);
    std::string pmat = "../../data/mouse54_grm.rel";
    std::string pmatid = "../../data/mouse54_grm.rel.id";
    EigenMat cor(pmat.c_str(), pmatid.c_str());
    for (int i = 0; i < cor.size1; ++i)
    {
        std::cout<<cor.rownames.eltOf(i)<< " ";
    }
    std::cout<<std::endl;
    cor.eigen();

    std::cout << "eigenvectors: "<<std::endl;
    for (int i = 0; i < cor.size1; ++i)
    {
        std::cout<<gsl_vector_get(cor.eigenvalues, i)<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"eigenvectors size1: "<<cor.eigenvectors->size1<<std::endl;
    std::cout<<"eigenvectors size2: "<<cor.eigenvectors->size2<<std::endl;
    std::cout << "eigenvectors: "<<std::endl;
    for (int i = 0; i < cor.eigenvectors->size1; ++i){
        for (int j=0; j < cor.eigenvectors->size2; ++j) {
            std::cout<< gsl_matrix_get(cor.eigenvectors, i, j) << " ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
    EigenMat * p = new EigenMat(pmat.c_str(), pmatid.c_str());
    std::cout<<"pointer test: "<<p<<std::endl;
    delete p;
    std::cout<<"delete pointer test: "<<p<<std::endl;
    p = &cor;
    std::cout<<"reassign pointer test: "<<p<<std::endl;
    p = nullptr;
    std::cout<<"set to NULL pointer test: "<<p<<std::endl;

}

TEST(MANOVA_READ, DISABLED_manova_test) {

    std::string pmat = "../../data/mouse54_grm.rel";
    std::string pmatid = "../../data/mouse54_grm.rel.id";
    std::string pstrain = "../../data/MPD_39504/trait.39504.txt";

    Dynum<string> strains;
    // read row names file
    ColumnReader rdr(pstrain.c_str(), (char *)"\t");
    int numtoks;
    while ((numtoks = rdr.getLine()) >= 0)
    {
        std::string strain_abbrev = rdr.getToken(0);
        int strIdx = strains.addElementIfNew(strain_abbrev);
        if (strIdx < 0)
        {
            std::cout << "Undefined strain abbrev: " << strain_abbrev << std::endl;
        }
    }
    //EXPECT_TRUE(aov);
    double P = 0,F=0;
    const char* pattern = "00000000000000000000000000101111";
    char* pat = strdup(pattern);
    int N = strlen(pattern);
    std::cout<<"pattern length: "<<N<<std::endl;
    makeUnprintable(pat);

    MANOVA aov(pmat.c_str(), pmatid.c_str(), 4);
    aov.setEigen();
    bool ok = aov.setNonQMarkMat(pat, strains);
    if (ok)
        aov.pillaiTrace(F, P);

    std::cout<<"FStat: "<<F<<std::endl;
    std::cout<<"Pvalue: "<<P<<std::endl;
    std::cout<<aov<<std::endl;
    EXPECT_GT(F,0);
    EXPECT_GT(P, 0);

}




