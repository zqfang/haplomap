#include <algorithm>
#include <functional>
#include "ghmap.h"
#include "gsl/gsl_cdf.h"


// defined globals
int numCategories = 1; // default for non-categorical data.
const int AACLASSES[] = {4, -1, 3, 1, 1, 4, 2, 0, 4, -1, 0, 4, 4, 2, -1, 4, 2, 0, 2, 2, -1, 4, 4, 5, 2, -1};
std::vector<std::string> catNames; // maps strIdx -> category name.
// Maps gene names to a string of A's, M's, and P's
std::unordered_map<std::string, std::string> geneExprMap(40000);
// Globals
std::unordered_map<std::string, GeneSummary *> geneTable; // for gene-oriented interface
std::vector<BlockSummary *> blocks; // global vector of all blocks.
int traceFStat = false;



// constructor
BlockSummary::BlockSummary(const char *chrnm, int num, int start, int size,
             int chrbeg, int chrend, const char *pat):
          blockIdx(num), blockStart(start), blockSize(size),
          chrBegin(chrbeg), chrEnd(chrend), isIgnored(false),
          FStat(INFINITY), pvalue(1.0), FDR(1.0), effect(0.0),
          relFStat(INFINITY), relPvalue(1.0), relFDR(1.0), relReject(false),
          numHaplo(-1), numInteresting(-1)
{
    chrName = strdup(chrnm);
    pattern = strdup(pat);
}

BlockSummary::~BlockSummary()
{
    /// FIXME: need to free memory if called strdup()
    free(pattern);
    free(chrName);
}

// summary of a block, read for the file.
std::string BlockSummary::updateCodonScore(std::string str)
  {
    int count = 0;
    size_t pos = 0;
    size_t endpos = 0;
    while (true)
    {
      pos = str.find("<", pos + 1);
      if (pos == std::string::npos)
        break;
      endpos = str.find(">", pos + 1);
      int aa1 = str[pos - 1] - 'A';
      int aa2 = str[endpos + 5] - 'A';
      if (AACLASSES[aa1] != AACLASSES[aa2])
      {
        count++;
        int X = 'X' - 'A';
        if (aa1 == X || aa2 == X)
        { // stop codon
          return "stop_codon";
        }
      }
    }
    if (str.find("SPLICE_SITE") != std::string::npos)
    {
      return "splicing";
    }
    if (count > 0)
    {
      return "codon_change";
    }
    return "no_codon_change";
  }

void updateGeneIsInteresting(BlockSummary *pb)
  { // BY
    for (std::map<std::string, std::string>::iterator giit = pb->geneIsCodingMap.begin(); giit != pb->geneIsCodingMap.end(); giit++)
    {
      if (giit->second == "0" || giit->second == "1")
      {
        pb->geneIsInteresting[giit->first] = "no_codon_change";
        continue;
      }
      pb->geneIsInteresting[giit->first] = pb->updateCodonScore(giit->second);
    }
  }

  // print a line of the blocks file.
  // blockIdx	blockStart	blockSize	chromosome	begin	end	pattern	pval	effect	genename genehascoding ...
void showBlockSum(std::ostream &os, bool isCategorical, BlockSummary *pb, std::vector<int> &strOrderVec)
  {
    os << pb->blockIdx;
    if (pb->isIgnored)
    {
      os << "\t(IGNORED)";
    }
    os << "\t" << pb->blockStart << "\t" << pb->blockSize
       << "\t" << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t";

    writeSortedPattern(os, pb->pattern, strOrderVec);

    os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) <<"\t" << pb->effect << "\t" <<pb->FDR;
    os << "\t" <<pb->relPvalue<<"\t"<<pb->relFDR<<pb->relReject;
    updateGeneIsInteresting(pb);

    // gene names and coding bits
    std::map<std::string, std::string> &geneIsCodingMap = pb->geneIsCodingMap;
    for (std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
    { ///EDITED BY
      os << "\t" << (*giit).first << "\t" << (*giit).second << "\t" << pb->geneIsInteresting[(*giit).first];
    }
    os << std::endl;
  }

// BY addition, output as such:
// gene	codon_flag	pattern	pval	effect	chromosome	begin	end	blockIdx	blockStart	blockSize	expression
void showGeneBlockByBlock(std::ostream &os, bool isCategorical, BlockSummary *pb, std::vector<int> &strOrderVec)
  {
    updateGeneIsInteresting(pb);
    std::map<std::string, std::string> &geneIsCodingMap = pb->geneIsCodingMap;
    std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin();
    if (giit == geneIsCodingMap.end())
    { // no genes in this block
      os << "None\tNone\t";
      writeSortedPattern(os, pb->pattern, strOrderVec);
      os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) <<"\t"<< pb->effect
         << "\t" << pb->FDR << "\t" << pb->relPvalue<<"\t"<<pb->relFDR<<"\t"<<pb->relReject
         << "\t" << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t"
         << pb->blockIdx << pb->blockStart << "\t" << pb->blockSize << std::endl;
    }
    else
    {
      for (std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
      {
        os << (*giit).first << "\t" << pb->geneIsInteresting[(*giit).first] << "\t";
        writeSortedPattern(os, pb->pattern, strOrderVec);
        os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) <<"\t"<< pb->effect
           << "\t" <<pb->FDR << "\t" << pb->relPvalue<<"\t"<<pb->relFDR<<"\t"<<pb->relReject
           << "\t" <<pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd
           << "\t" <<pb->blockIdx << pb->blockStart << "\t" << pb->blockSize;
        std::string gname = (*giit).first;
        upcase(gname);
        if (geneExprMap.find(gname) == geneExprMap.end())
        {
          os << "\t-------------" << std::endl;
        }
        else
        {
          os << "\t" << geneExprMap[gname] << std::endl;
        }
      }
    }
  }



void showBlockSums(std::ostream &os, bool isCategorical,
        std::vector<BlockSummary *> &blocks, float cutoff, std::vector<int> &strOrderVec)
{
    for (unsigned i = 0; i < blocks.size(); i++)
    {
      if (!isCutoff(isCategorical, cutoff, blocks[i]))
      {
        showBlockSum(os, isCategorical, blocks[i], strOrderVec);
      }
    }
}

void showGeneBlockByBlocks(std::ostream &os, bool isCategorical,
        std::vector<BlockSummary *> &blocks, float cutoff, std::vector<int> &strOrderVec)
  {
    for (unsigned i = 0; i < blocks.size(); i++)
    {
      if (!isCutoff(isCategorical, cutoff, blocks[i]))
      {
        showGeneBlockByBlock(os, isCategorical, blocks[i], strOrderVec);
      }
    }
  }


bool BlocksComparator::operator()(const BlockSummary *pb1, const BlockSummary *pb2) const
  {
    // Use F statistic or p value, depending on whether it's categorical or not.
    if (isCategorical)
    {
      if (pb1->FStat > pb2->FStat)
      {
        return true;
      }
      else if (pb1->FStat < pb2->FStat)
      {
        return false;
      }
    }
    else
    {
      if (pb1->pvalue < pb2->pvalue)
      {
        return true;
      }
      else if (pb1->pvalue > pb2->pvalue)
      {
        return false;
      }
    }

    // Break ties by choosing blocks with the fewest haplotypes.
    int nh1 = pb1->numHaplo; //numHaplotypes(pb1->pattern);
    int nh2 = pb2->numHaplo; //numHaplotypes(pb2->pattern);

    if (nh1 < nh2)
    {
      return true;
    }
    else if (nh2 < nh1)
    {
      return false;
    }

    // Then order by SNP changes
    int snpChangeScore1 = pb1->numInteresting; //interestingChanges(pb1->geneIsCodingMap);
    int snpChangeScore2 = pb2->numInteresting; //interestingChanges(pb2->geneIsCodingMap);

    if (snpChangeScore1 < snpChangeScore2)
    {
      return false;
    }
    else if (snpChangeScore2 < snpChangeScore1)
    {
      return true;
    }
    // Then order by chromosome
    char *chr1 = pb1->chrName;
    char *chr2 = pb2->chrName;

    if (strcmp(chr1, chr2) != 0)
    {
      // deal with complexities of comparing chr names.
      // put "M" last
      if (strcmp(chr2, "M") == 0)
      {
        return true;
      }
      else if (strcmp(chr1, "M") == 0)
      {
        return false;
      }
      // Y next to last
      else if (strcmp(chr2, "Y") == 0)
      {
        return true;
      }
      else if (strcmp(chr1, "Y") == 0)
      {
        return false;
      }
      // X just before Y
      else if (strcmp(chr2, "X") == 0)
      {
        return true;
      }
      else if (strcmp(chr1, "X") == 0)
      {
        return false;
      }

      // otherwise, in numerical order
      int numChr1 = std::atoi(chr1); // char* to integer
      int numChr2 = std::atoi(chr2);
      if (numChr1 < numChr2)
      {
        return true;
      }
      else if (numChr1 > numChr2)
      {
        return false;
      }
      else
      {
        return false;
      }
    }
    // chromosome names were same, so use position
    else if (pb1->chrBegin < pb2->chrBegin)
    {
      return true;
    }
    else
    {
      return false;
    }
  };
void readBlockSummary(char *fname, char *geneName, bool ignoreDefault)
{
    ColumnReader rdr(fname, (char *)"\t");

    int numtoks;
    while ((numtoks = rdr.getLine()) >= 0)
    {
        // file has "Chromosome\tBlockNum\tBlockstartSNPnum\tSize\tChrbegin\tChrend\tPattern\tgene1\tcodon1\t...\n"

        // If geneName is defined, only process the lines containing that gene name in the gene name list,
        // unless geneName = "*", in which case it reads all blocks.

        if (!geneName || (std::find(rdr.begin() + 7, rdr.end(), geneName) != rdr.end()))
        {

            BlockSummary *pBlock = new BlockSummary(rdr.getToken(0).c_str(),
                                                    std::stoi(rdr.getToken(1)),
                                                    std::stoi(rdr.getToken(2)),
                                                    std::stoi(rdr.getToken(3)),
                                                    std::stoi(rdr.getToken(4)),
                                                    std::stoi(rdr.getToken(5)),
                                                    rdr.getToken(6).c_str());
            // Add to blocks.
            blocks.push_back(pBlock);
            BlockSummary *pLastBlock = blocks.back();
            // For each gene, add gene summary to the gene table (if not there already),
            // and insert block in the block set.
            // This iterates over list of alternating gene name/codon tokens.
            for (std::vector<std::string>::iterator git = rdr.begin() + 7; git != rdr.end(); git += 2)
            {
                // *(git+1) is codon for gene name (*git).  Convert from string to bool using cond expr.
                pLastBlock->geneIsCodingMap[*git] = *(git + 1);
                std::unordered_map<std::string, GeneSummary *>::iterator entIt = geneTable.find(*git);

                if (entIt == geneTable.end())
                {
                    // create new entry.
                    // cout << "  Creating new gene summary for " << *git << endl;
                    geneTable[*git] = new GeneSummary(ignoreDefault);
                    geneTable[*git]->name = *git;
                }
                // Add block to GeneSummary
                //      cout << "  Adding block for " << *git << ": " << pLastBlock->blockIdx << endl;
                geneTable[*git]->blocks.push_back(pLastBlock);
            }

            // cout << "Gene names for block " << blocks.back().blockIdx << ": " << blocks.back().geneNames << endl;
        }
    }

    // Separate pass to fix up patterns.
    // patterns are all the same length.
    // FIXME: this causes a seg fault when there are no blocks (happens under funny conditions).
    numStrains = strlen(blocks[0]->pattern);

    for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
    {
        BlockSummary *block = blocks[blkIdx];
        makeUnprintable(block->pattern);
    }
}

void showGeneBlockSums(std::ostream &os, bool isCategorical, std::vector<BlockSummary *> &blocks,
                       float cutoff, std::vector<int> &strOrderVec, std::vector<GeneSummary *> genesList)
{
    std::vector<std::string> genesOver; // all genes that have been already printed
    for (unsigned i = 0; i < blocks.size(); i++)
    {
        if (!isCutoff(isCategorical, cutoff, blocks[i]))
        {
            BlockSummary *pb = blocks[i];
            std::map<std::string, std::string> &geneIsCodingMap = pb->geneIsCodingMap;
            std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin();
            if (geneIsCodingMap.size() == 0)
            { // no genes in this block
                updateGeneIsInteresting(pb);
                os << "None\tNone\t";
                writeSortedPattern(os, pb->pattern, strOrderVec);
                os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) << "\t" << pb->effect << "\t";
                os << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t";
                os << pb->blockIdx << pb->blockStart << "\t" << pb->blockSize << std::endl;
            }
            else
            {
                for (; giit != geneIsCodingMap.end(); giit++)
                {
                    std::string gname = giit->first;
                    if (std::find(genesOver.begin(), genesOver.end(), gname) != genesOver.end())
                    {
                        continue;
                    }
                    genesOver.push_back(gname);
                    for (std::vector<GeneSummary *>::iterator git = genesList.begin(); git != genesList.end(); git++)
                    {
                        if ((*git)->name == gname)
                        {
                            for (unsigned j = 0; j < ((*git)->blocks).size(); j++)
                            {
                                BlockSummary *pb_t = (*git)->blocks[j];
                                updateGeneIsInteresting(pb_t);
                                //os << gname << "\t" << pb_t->geneIsInteresting[gname] << "\t";
                                os << gname << "\t" << pb_t->geneIsCodingMap[gname] << "\t";
                                writeSortedPattern(os, pb_t->pattern, strOrderVec);
                                os << "\t" << (isCategorical ? pb_t->FStat : pb_t->pvalue) <<"\t"<< pb_t->effect;
                                os << "\t" << pb_t->FDR;
                                os << "\t" << pb_t->relPvalue<<"\t"<<pb_t->relFDR <<"\t"<<pb_t->relReject<<"\t";
                                os << pb_t->chrName << "\t" << pb_t->chrBegin << "\t" << pb_t->chrEnd << "\t";
                                os << pb_t->blockIdx << pb_t->blockStart << "\t" << pb_t->blockSize;
                                std::string upname = gname;
                                upcase(upname);
                                if (geneExprMap.find(upname) == geneExprMap.end())
                                {
                                    os << "\t-------------" << std::endl;
                                }
                                else
                                {
                                    os << "\t" << geneExprMap[upname] << std::endl;
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
}

void writeGeneBlockSums(bool isCategorical, char *outputFileName, char *datasetName,
                        std::vector<std::vector<float>> &phenvec, std::vector<BlockSummary *> &blocks, float pvalueCutoff)
{
    std::ofstream blockout(outputFileName);
    if (!blockout.is_open())
    {
        std::cout << "Open of file \"" << outputFileName << "\" failed: ";
        perror("");
        exit(1);
    }
    // Datasetname
    blockout << datasetName << std::endl;
    // sort strains by phenvec value
    std::vector<int> strOrderVec(numStrains); // will contain strain indices.
    sortStrainsByPheno(phenvec, strOrderVec);

    // output strain names.
    std::vector<int>::iterator stoEnd1 = strOrderVec.end();
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        blockout << strainAbbrevs.eltOf(str1);
        if (stoIt1 + 1 < stoEnd1)
        {
            blockout << "\t";
        }
    }
    blockout << std::endl;

    // output phenotype values.
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        if (phenvec[str1].size() > 1) {
            for (unsigned i = 0; i < phenvec[str1].size() - 1; i++)
                blockout << phenvec[str1][i] << ",";
        }
        blockout << phenvec[str1].back();

        if (stoIt1 + 1 < stoEnd1)
        {
            blockout << "\t";
        }
    }
    blockout << std::endl;

    GenesComparator gcomp(isCategorical);
    std::vector<GeneSummary *> genes;
    genes.reserve(geneTable.size());
    //transform(geneTable.begin(), geneTable.end(), back_inserter(genes), select2nd<std::unordered_map<string, GeneSummary *> >());
    transform(geneTable.begin(), geneTable.end(), back_inserter(genes),
              std::bind(&std::unordered_map<std::string, GeneSummary *>::value_type::second, std::placeholders::_1 ));

    sort(genes.begin(), genes.end(), gcomp);
    // write header
//    blockout << "gene_name\tcondon\tpattern\t";
//    blockout << (isCategorical ? "FStat" : "pvalue");
//    blockout << "\teffect\tFDR\tpopPvalue\tpopFDR\tpop\tchr\tstart\tend\t";
//    blockout << "blockIdx\tblockSize\tgene_expression\n";
    showGeneBlockSums(blockout, isCategorical, blocks, pvalueCutoff, strOrderVec, genes);
}


void writeGeneBlockByBlocks(bool isCategorical, char *outputFileName, char *datasetName,
                            std::vector<std::vector<float>> &phenvec,std::vector<BlockSummary *> &blocks, float pvalueCutoff)
{
    std::ofstream blockout(outputFileName);
    if (!blockout.is_open())
    {
        std::cout << "Open of file \"" << outputFileName << "\" failed: ";
        perror("");
        exit(1);
    }
    // Datasetname
    blockout << datasetName << std::endl;
    // sort strains by phenvec value
    std::vector<int> strOrderVec(numStrains); // will contain strain indices.
    sortStrainsByPheno(phenvec, strOrderVec);

    // output strain names.
    std::vector<int>::iterator stoEnd1 = strOrderVec.end();
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        blockout << strainAbbrevs.eltOf(str1);
        if (stoIt1 + 1 < stoEnd1)
        {
            blockout << "\t";
        }
    }
    blockout << std::endl;

    // output phenotype values.
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        if (phenvec[str1].size() > 1) {
            for (unsigned i = 0; i < phenvec[str1].size() - 1; i++)
                blockout << phenvec[str1][i] << ",";
        }
        blockout << phenvec[str1].back();
        if (stoIt1 + 1 < stoEnd1)
        {
            blockout << "\t";
        }
    }
    blockout << std::endl;

    // blockout << "gene_name\tcondon\tpattern\t";
    // blockout << (isCategorical ? "FStat" : "pvalue");
    // blockout << "\teffect\tFDR\tpopPvalue\tpopFDR\tpop\tchr\tstart\tend\t";
    // blockout << "blockIdx\tblockSize\tgene_expression\n";
    showGeneBlockByBlocks(blockout, isCategorical, blocks, pvalueCutoff, strOrderVec);
}

void writeBlockSums(bool isCategorical, char *outputFileName,
                    char *datasetName, std::vector<std::vector<float> > &phenvec,
                    std::vector<BlockSummary *> &blocks, float pvalueCutoff)
{
    std::ofstream blockout(outputFileName);
    if (!blockout.is_open())
    {
        std::cout << "Open of file \"" << outputFileName << "\" failed: ";
        perror("");
        exit(1);
    }

    // Datasetname
    blockout << datasetName << std::endl;

    // sort strains by phenvec value

    std::vector<int> strOrderVec(numStrains); // will contain strain indices.
    sortStrainsByPheno(phenvec, strOrderVec);

    // output strain names.
    // FIXME:  Start using vecfuns written for Ravi microarray analysis
    std::vector<int>::iterator stoEnd1 = strOrderVec.end();
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        blockout << strainAbbrevs.eltOf(str1);
        if (stoIt1 + 1 < stoEnd1)
        {
            blockout << "\t";
        }
    }
    blockout << std::endl;

    // output phenotype values.
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        if (isCategorical)
        {
            blockout << catNames[str1];
        }
        else
        {
            if (phenvec[str1].size() > 1) {
                for (unsigned i = 0; i < phenvec[str1].size() - 1; i++)
                    blockout << phenvec[str1][i] << ",";
            }
            blockout << phenvec[str1].back();
        }
        if (stoIt1 + 1 < stoEnd1)
        {
            blockout << "\t";
        }
    }
    blockout << std::endl;
//    // write header
//    blockout <<"blockIdx\tblockStart\tblockSize\tchr\tstart\tend\tpattern\t";
//    blockout <<(isCategorical ? "FStat":"pvalue");
//    blockout <<"\teffect\tFDR\tpopPvalue\tpopFDR\tpop\tgene_name\tmutate\tcondon\n";

    showBlockSums(blockout, isCategorical, blocks, pvalueCutoff, strOrderVec);
}

// Write gene-oriented summary.
void writeGeneSums(bool isCategorical, char *outputFileName,
                   char *datasetName, std::vector<std::vector<float> > &phenvec,
                   std::vector<BlockSummary *> &blocks, float cutoff, bool filterCoding)
{
    std::ofstream genesout(outputFileName);
    if (!genesout.is_open())
    {
        std::cout << "Open of file \"" << outputFileName << "\" failed: ";
        perror("");
        exit(1);
    }

    // Datasetname
    genesout << datasetName << std::endl;

    // sort strains by phenvec value

    std::vector<int> strOrderVec(numStrains); // will contain strain indices.
    sortStrainsByPheno(phenvec, strOrderVec);

    // output strain names.
    // FIXME:  Start using vecfuns written for Ravi microarray analysis
    //  ... or STL algorithms!
    std::vector<int>::iterator stoEnd1 = strOrderVec.end();
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        genesout << strainAbbrevs.eltOf(str1);
        if (stoIt1 + 1 < stoEnd1)
        {
            genesout << "\t";
        }
    }
    genesout << std::endl;

    // output phenotype values.
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        if (isCategorical)
        {
            genesout << catNames[str1];
        }
        else
        {
            if (phenvec[str1].size() > 1) {
                for (unsigned i = 0; i < phenvec[str1].size() - 1; i++)
                    genesout << phenvec[str1][i] << ",";
            }
            genesout << phenvec[str1].back();
        }
        if (stoIt1 + 1 < stoEnd1)
        {
            genesout << "\t";
        }
    }
    genesout << std::endl;

    GenesComparator gcomp(isCategorical);

    // Copy genesTable values into a vector and sort using GenesComparator
    std::vector<GeneSummary *> genes;
    genes.reserve(geneTable.size());
    std::transform(geneTable.begin(), geneTable.end(), back_inserter(genes),
              std::bind(&std::unordered_map<std::string, GeneSummary *>::value_type::second, std::placeholders::_1 ));
    std::sort(genes.begin(), genes.end(), gcomp);
    // write header
//    genesout <<"gene_name\tcondon\tpattern\t";
//    genesout << (isCategorical ? "FStat" : "pvalue");
//    genesout << "\teffect\tFDR\tpopPvalue\tpopFDR\tpop\tchr\tstart\tend\tgene_expression\n";
    // write them out.
    for (std::vector<GeneSummary *>::iterator git = genes.begin(); git != genes.end(); git++)
    {

        if ((*git)->isIgnored)
            continue; //simply don't write these genes
        BlockSummary *pBestBlock = (*git)->blocks[0];
        if (isCategorical)
        {
            if (pBestBlock->FStat < cutoff)
            {
                break;
            }
        }
        else
        {
            if (pBestBlock->pvalue > cutoff)
            {
                break;
            }
        }

        // write gene name and its coding bit.
        std::string &gname = (*git)->name;
        std::string ugname = gname; // gene names in expression data are upper case.
        upcase(ugname);

        // coding bit is subtle.  This iterates over all of the blocks that are < pvalue cutoff
        // (or > for FStat).  If any has a non-synonmous SNP, the bit is 1.  So, the coding
        // change may not be in the best block for the gene, and may change if the pvalue cutoff
        // is changed.
        //string isCodingStr = "";
        bool isCoding = false;
        bool hasInteresting = false;
        bool hasSpliceChange = false;
        //ofstream debug_log;
        //debug_log.open("debug.log",ios::app);
        for (std::vector<BlockSummary *>::iterator blit = (*git)->blocks.begin(); blit != (*git)->blocks.end(); blit++)
        {
            if (isCutoff(isCategorical, cutoff, *blit))
            {
                break;
            }
            else
            {
                //debug_log << gname << "\t" << (*blit)->geneIsCodingMap[gname] << endl;
                if ((*blit)->geneIsCodingMap[gname] != "0")
                { //if the gene had SNPs marked as NON_SYNONYMOUS_CODING, with <->, or as SPLICE_SITE, isCoding is true
                    isCoding = true;
                }
                hasInteresting |= ((*blit)->numInteresting > 0);                                            //has a major amino acid change
                hasSpliceChange |= (((*blit)->geneIsCodingMap[gname]).find("SPLICE_SITE") != std::string::npos);
                //if "SPLICE_SITE" was in there that means that the gene had a splice change
            }
        }
        //New thing: -1 means not coding, 0 means coding but not important
        //Anything else is the number of important
        int codingCode = -1;
        if (isCoding)
            codingCode = 0;
        if (hasInteresting)
            codingCode = 1;
        if (hasSpliceChange)
            codingCode = 2; //if this is true, will overwrite codingCode=1
        updateGeneIsInteresting(pBestBlock);
        if ((isCoding || !filterCoding) || hasSpliceChange)
        {
            //genesout << gname << "\t" << pBestBlock->geneIsInteresting[(*git)->name] << "\t";
            //genesout << gname << "\t" << pBestBlock->geneIsCodingMap[gname] << "\t";
            genesout << gname << "\t" << codingCode << "\t";
            writeSortedPattern(genesout, pBestBlock->pattern, strOrderVec);
            genesout << "\t" << (isCategorical ? pBestBlock->FStat : pBestBlock->pvalue);
            genesout << "\t" << pBestBlock->effect <<"\t"<< pBestBlock->FDR;
            genesout << "\t" << pBestBlock->relPvalue<<"\t"<< pBestBlock->relFDR << "\t"<<pBestBlock->relReject;
            genesout << "\t" << pBestBlock->chrName << "\t" << pBestBlock->chrBegin << "\t" << pBestBlock->chrEnd;

            // write gene expression values
            if (geneExprMap.find(ugname) == geneExprMap.end())
            {
                // no data for that gene name
                genesout << "\t-------------" << std::endl;
            }
            else
            {
                genesout << "\t" << geneExprMap[ugname] << std::endl;
            }
        }
    }
}


void sortStrainsByPheno(std::vector<std::vector<float>> &phenvec, std::vector<int> &strOrderVec)
{
    std::vector<int>::iterator stoEnd = strOrderVec.end();
    int i = 0;
    for (std::vector<int>::iterator stoIt = strOrderVec.begin(); stoIt != stoEnd; stoIt++)
    {
        *stoIt = i++;
    }
    // This will sort lexicographically, which is the right thing.
    // For categorical values, we just want equal values together.
    IndexComparator<std::vector<float>, std::less<std::vector<float>>> idxCompare(&phenvec);
    std::stable_sort(strOrderVec.begin(), strOrderVec.end(), idxCompare);
}

// Count the number of defined strains
int numDefinedStrains(char *pattern)
{
    int count = 0;
    for (int str1 = 0; str1 < numStrains; str1++)
    {
        if (pattern[str1] != '?')
        {
            count++;
        }
    }
    return count;
}

// returns maximum eq. class + 1.
int numHaplotypes(char *pattern)
{
    int numHap = -1;
    for (int str1 = 0; str1 < numStrains; str1++)
    {
        int hap = pattern[str1];
        if (pattern[str1] != '?' && numHap < hap)
        {
            numHap = hap;
        }
    }
    return numHap + 1;
}

/* Returns a score that represents the interestingness of codon changes*/
int scoreChanges(std::string str)
{
    int count = 0;
    size_t pos = 0;
    size_t endpos = 0;
    while (true)
    {
        pos = str.find("<", pos + 1);
        if (pos == std::string::npos)
            break;
        endpos = str.find(">", pos + 1);
        int aa1 = str[pos - 1] - 'A';
        int aa2 = str[endpos + 5] - 'A';
        if (AACLASSES[aa1] != AACLASSES[aa2])
        {
            count++;
            int X = 'X' - 'A';
            if (aa1 == X || aa2 == X)
                count += 2; // Higher weight if stop codon
        }
    }
    return count;
}

int interestingChanges(const std::map<std::string, std::string> &geneCodingMap)
{
    int changeCount = 0;
    for (std::map<std::string, std::string>::const_iterator git = geneCodingMap.begin();
         git != geneCodingMap.end(); git++)
    {
        if (git->second == "0" || git->second == "1")
            continue;
        changeCount += scoreChanges(git->second);
//        changeCount+=countInStr(git->second, "<->");
//        changeCount+=countInStr(git->second, "SPLICE_SITE");
    }
    return changeCount;
}

// Mark each block as ignored unless it has a coding gene.
void filterCodingBlocks()
{
    for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
    {
        BlockSummary *pBlock = blocks[blkIdx];
        // look for a coding gene
        bool keep = false;
        std::map<std::string, std::string> &gicMap = pBlock->geneIsCodingMap;

        for (std::map<std::string, std::string>::iterator gicmit = gicMap.begin(); gicmit != gicMap.end(); gicmit++)
        {
            if ((*gicmit).second != "0")
            {
                keep = true;
                break;
            }
        }
        pBlock->isIgnored |= !keep; // |= as style to allow multiple filters.
    }
}

void filterEqualBlocks(std::vector<int> equalRegions)
{
    for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
    {
        BlockSummary *pBlock = blocks[blkIdx];
        // look for a coding gene
        bool keep = true;
        char allele = -1;
        for (unsigned int i = 0; i < equalRegions.size(); i++)
        {
            if (allele == -1)
            {
                allele = pBlock->pattern[equalRegions[i]];
                continue;
            }
            keep &= (allele == (pBlock->pattern)[equalRegions[i]]);
        }
        pBlock->isIgnored |= !keep; // |= as style to allow multiple filters.
    }
}

// Read the file of gene expression data
void readCompactGeneExpr(char *fname)
{
    ColumnReader rdr(fname, (char *)"\t");
    int numtoks;
    while ((numtoks = rdr.getLine()) >= 0)
    {
        // A typical line: "Myc	PAPAAAAAAAPAA"
        std::string geneName = rdr.getToken(0);
        std::string present = rdr.getToken(1);

        upcase(geneName);

        // Gene names in expr file are all upper case, so use upper case comparisons.
        // Don't really need to do upcase here, but doing it anyway in case someone
        // substitutes a file with mixed-case names.
        geneExprMap[geneName] = present;
    }
}

// Read a file of quantitative phenotypes.
void readQPhenotypes(char *fname, std::vector<std::vector<float>> &phenvec)
{
    ColumnReader rdr(fname, (char *)"\t");

    int numtoks;

    while ((numtoks = rdr.getLine()) >= 0)
    {
        // file has "Abbrev\tValue\n"
        if (numtoks != 2)
        {
            std::cout << "Warning: numtoks = " << numtoks << std::endl;
        }

        // FIXME: some unnecessary string copies
        std::string strain_abbrev = rdr.getToken(0);
        std::vector<float> qphen;
        //qphen.push_back(std::stof(rdr.getToken(1)));
        //    int strIdx = strainAbbrevs.hasIndex(strain_abbrev);
        int strIdx = strainAbbrevs.addElementIfNew(strain_abbrev);
        if (strIdx < 0)
        {
           std::cout << "Undefined strain abbrev: " << strain_abbrev << std::endl;
        }
        //phenvec[strIdx] = qphen;
        // MARK: handle same animal with multiple values
        phenvec[strIdx].push_back(std::stof(rdr.getToken(1)));
    }
}

void setBlockStats()
{
    for (unsigned blkIdx = 0; blkIdx < blocks.size(); blkIdx++)
    {
        BlockSummary *pBlock = blocks[blkIdx];
        //if(!pBlock->isIgnored) {
        pBlock->numHaplo = numHaplotypes(pBlock->pattern);
        pBlock->numInteresting = interestingChanges(pBlock->geneIsCodingMap);
    }
}
// Read a file of categorical phenotypes.
void readCPhenotypes(char *fname, std::vector<std::vector<float>> &phenvec)
{

    ColumnReader rdr(fname, (char *)"\t");
    Dynum<std::string> categories; // assigned the distinct categories consecutive indices, starting at 0.

    catNames.resize(numStrains);

    int numtoks;

    while ((numtoks = rdr.getLine()) >= 0)
    {
        // file has "Abbrev\t\Value"
        if (numtoks != 2)
        {
            std::cout << "Warning: numtoks = " << numtoks << std::endl;
        }

        // FIXME: some unnecessary string copies
        std::string strain_abbrev = rdr.getToken(0);
        std::string catname = rdr.getToken(1);
        int strIdx = strainAbbrevs.addElementIfNew(strain_abbrev);
        catNames[strIdx] = catname;
        categories.addElementIfNew(catname);
    }

    // build category vectors.
    numCategories = categories.size();

    for (int strIdx = 0; strIdx < numStrains; strIdx++)
    {
        // build vector value for this category and store in phenvec.
        std::vector<float> cphen(categories.size(), 0.0F); // initialize to 0.
        cphen[categories.indexOf(catNames[strIdx])] = 1.0F;
        phenvec[strIdx] = cphen;
    }

    if (traceFStat)
    {
        std::cout << "Phenotype values" << std::endl;
        categories.dump();

        std::cout << "Phenotype vectors for strains" << std::endl;
        for (int str = 0; str < numStrains; str++)
        {

            for(auto &t: phenvec[str])
                std::cout << t << " ";
        }
        std::cout << std::endl;
    }
}

//read in equal class
void readEqualFile(char *fname, std::vector<int> &equalStrains)
{
    ColumnReader rdr(fname, (char *)"\t");
    //numStrains

    int numtoks;
    int i = 0;
    while ((numtoks = rdr.getLine()) >= 0)
    {
        std::string eqclass = rdr.getToken(1);
        if (eqclass != "0")
            equalStrains.push_back(i);
        i++;
    }
}

//reads first column of tsv fname into vector vec
void readFileToVec(char *fname, std::vector<std::string> &vec)
{
    ColumnReader rdr(fname, (char *)"\t");
    //numStrains

    int numtoks;
    while ((numtoks = rdr.getLine()) >= 0)
    {
        vec.push_back(rdr.getToken(0));
    }
}

void filterGoTerms(char *fname, std::vector<std::string> terms)
{
    ColumnReader rdr(fname, (char *)"\t");
    int numtoks;
    std::vector<std::string>::iterator startT = terms.begin();
    std::vector<std::string>::iterator endT = terms.end();
    while ((numtoks = rdr.getLine()) >= 0)
    {
        std::string geneName = rdr.getToken(0);
        std::unordered_map<std::string, GeneSummary *>::iterator it = geneTable.find(geneName);
        if (it == geneTable.end())
        {
            continue;
        }
        for (int i = 1; i < numtoks; i++)
        {
            if (std::find(startT, endT, rdr.getToken(i)) != endT)
            {
                geneTable[geneName]->isIgnored = false;
            }
            geneTable[geneName]->goTerms.push_back(rdr.getToken(i));
        }
    }
}

// renumber eqclasses in pattern so the increase from left to right
void writeSortedPattern(std::ostream &os, char *pattern, std::vector<int> &strOrderVec)
{
    int numHaplo = numHaplotypes(pattern);
    char *sortedEqMap = (char *)malloc(numHaplo);
    memset(sortedEqMap, '?', numHaplo);
    char curEq = 0;

    for (int strIdx = 0; strIdx < numStrains; strIdx++)
    {
        int str = strOrderVec[strIdx];
        int eqclass = (int)pattern[str];
        if ('?' != eqclass)
        {
            if ('?' == sortedEqMap[eqclass])
            {                                 // class hasn't been mapped yet.
                sortedEqMap[eqclass] = curEq++; // renumber for left-to-right color consistency
            }
            os << (char)(sortedEqMap[eqclass] + '0');
        }
        else
        {
            os << '?';
        }
    }
    // memory leak?
    free(sortedEqMap);
}

// Some vector arithmetic.
// destroys first argument (like +=)
void addVectors(std::vector<float> &v1, std::vector<float> &v2)
{
    if (v1.size() != v2.size())
    {
        std::cout << "addVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() << std::endl;
        exit(1);
    }
    std::vector<float>::iterator vend = v1.end();
    std::vector<float>::iterator vit2 = v2.begin();
    for (std::vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
    {
        *vit1 += *vit2;
        vit2++;
    }
}

void subtractVectors(std::vector<float> &v1, std::vector<float> &v2)
{
    if (v1.size() != v2.size())
    {
        std::cout << "subtractVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() <<std::endl;
        exit(1);
    }
    std::vector<float>::iterator vend = v1.end();
    std::vector<float>::iterator vit2 = v2.begin();
    for (std::vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
    {
        *vit1 -= *vit2;
        vit2++;
    }
}

//
float dotVectors(std::vector<float> &v1, std::vector<float> &v2)
{
    float result = 0.0;
    if (v1.size() != v2.size())
    {
        std::cout << "dotVectors:  Vector sizes differ: " << v1.size() << " vs. " << v2.size() << std::endl;
        exit(1);
    }
    std::vector<float>::iterator vend = v1.end();
    std::vector<float>::iterator vit2 = v2.begin();
    for (std::vector<float>::iterator vit1 = v1.begin(); vit1 < vend; vit1++)
    {
        result += (*vit1) * (*vit2);
        vit2++;
    }
    return result;
}

// multiply by scalar.  Destroys first argument.
void scaleVector(std::vector<float> &v1, float c)
{
    std::vector<float>::iterator vend = v1.end();
    for (std::vector<float>::iterator vit = v1.begin(); vit < vend; vit++)
    {
        *vit *= c;
    }
}

// convert digits 0-9 in string to \000..\011 (but leave '?' printable).
/// Mark: make ascii starts from 0, not '0'.
void makeUnprintable(char *pattern)
{
    char *p = pattern;
    while (*p != 0)
    {
        if (*p != '?')
        {
            *p -= '0';
        }
        p++;
    }
}


void bh_fdr(std::vector<BlockSummary *> & pval, float alpha, bool flag)
{
    //std::vector<bool> reject(pval.size(), false);
    float m = pval.size();
    uint32_t k = pval.size(); // This is the rank, doesn't need to be double.
    float factor;
    float p;
    float previous_fdr;
    //BlockSummary* pBlock = new BlockSummary();
    // stored padj
    if (flag) {
        std::stable_sort(pval.begin(), pval.end(),
                         [](BlockSummary* x, BlockSummary* y) {return x->pvalue > y->pvalue;});
        previous_fdr =1.0;
        for (unsigned i = 0; i < pval.size(); ++i) {
            factor = k / m;
            p = pval[i]->pvalue;
            //if (p <= factor * alpha) {
            //    pval[i]->relReject = true;
            //}
            p /= factor;
            pval[i]->FDR = std::min(p, previous_fdr); // accumulate minimum
            previous_fdr = pval[i]->FDR;
            k--; //Decrease rank
        }
//        std::accumulate(pval.begin(), pval.end(), pBlock,
//                        [](BlockSummary* x, BlockSummary* y)
//                        { return std::min(x->FDR, y->FDR); });
    } else {
        std::stable_sort(pval.begin(), pval.end(),
                         [](BlockSummary* x, BlockSummary* y) {return x->relPvalue > y->relPvalue;});
        previous_fdr = 1.0;
        for (unsigned i = 0; i < pval.size(); ++i) {
            factor = k / m;
            p = pval[i]->relPvalue;
            if (p <= factor * alpha) {
                pval[i]->relReject = true;
            }
            p /= factor;
            pval[i]->relFDR = std::min(p, previous_fdr);
            previous_fdr = pval[i]->relFDR;
            k--; //Decrease rank
        }
//        std::accumulate(pval.begin(), pval.end(), pBlock->mFDR,
//                        [](BlockSummary* x, BlockSummary* y)
//                                   { return std::min(x->mFDR, y->mFDR) ; });

    }
    //delete pBlock;
}

