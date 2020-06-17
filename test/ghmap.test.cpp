//
// Created by Zhuoqing Fang on 6/16/20.
//

#include <string>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "ghmap.h"


using namespace std;
TEST(READ_PHENO, cat_and_quant_read) {
    int numStrains = 12;
    vector<vector<float>> phenvec(numStrains);
    bool isCategorical = false;
    char* path = (char *)"../../data/trait.0000.txt";
    if (isCategorical)
    {
        readCPhenotypes(path, phenvec);
    }
    else
    {
        readQPhenotypes(path, phenvec);
    }

    cout << "Phenotype vectors for strains" << endl;
    for (int str = 0; str < numStrains; str++)
    {  
        cout << phenvec[str][0] << " ";
    }
    cout << endl;

}

TEST(READ_BLOCK, DISABLED_read_block_summary) {
    char* path = (char *)"../../data/chrX.hblocks.txt";
    std::vector<BlockSummary *> blocks; 
    readBlockSummary(path, NULL, false);
    BlockSummary *pBlock = blocks.back();
    EXPECT_NE(pBlock, nullptr);
    cout << "Block: " << pBlock->chrName << "\t" << pBlock->blockIdx << "\t"<<endl;
    showPattern(pBlock->pattern);

}


