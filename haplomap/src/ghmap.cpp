#include <algorithm>
#include <functional>
#include "ghmap.h"

// defined globals
int numCategories = 1; // default for non-categorical data.
const int AACLASSES[] = {4, -1, 3, 1, 1, 4, 2, 0, 4, -1, 0, 4, 4, 2, -1, 4, 2, 0, 2, 2, -1, 4, 4, 5, 2, -1};
std::vector<std::string> catNames; // maps strIdx -> category name.
// Maps gene names to a string of A's, M's, and P's
std::unordered_map<std::string, std::string> geneExprMap(40000);
// Globals
std::unordered_map<std::string, GeneSummary *> geneTable; // for gene-oriented interface
std::vector<BlockSummary *> blocks; // global vector of all blocks.
std::vector<std::string> geneExprHeader;
int traceFStat = false;

// std::unordered_map<std::string, int> PRIOR = {{"HIGH", 2}, {"MODERATE", 1}, {"LOW", 0}, {"MODIFIER", -1}};
std::unordered_map<std::string, int> CSQs = { 
    {"transcript_ablation", 2},
    {"splice_acceptor_variant",2},
    {"splice_donor_variant", 2},
    {"stop_gained",2},
    {"frameshift_variant", 2},
    {"stop_lost", 2},
    {"start_lost", 2},
    {"transcript_amplification", 2},
    {"inframe_insertion", 1},
    {"inframe_deletion", 1},
    {"missense_variant", 1},
    {"protein_altering_variant", 1},
    {"splice_region_variant", 0},
    {"incomplete_terminal_codon_variant", 0},
    {"start_retained_variant", 0},
    {"stop_retained_variant", 0},
    {"synonymous_variant", 0},
    {"coding_sequence_variant", -1},
    {"mature_miRNA_variant", -1},
    {"5_prime_UTR_variant", -1},
    {"3_prime_UTR_variant",-1},
    {"non_coding_transcript_exon_variant",-1},
    {"intron_variant", -1},
    {"NMD_transcript_variant", -1},
    {"non_coding_transcript_variant", -1},
    {"upstream_gene_variant", -1},
    {"downstream_gene_variant", -1},
    {"TFBS_ablation", -1},
    {"TFBS_amplification", -1},
    {"TF_binding_site_variant", -1},
    {"regulatory_region_ablation", -1},
    {"regulatory_region_amplification", -1},
    {"feature_elongation", -1},
    {"regulatory_region_variant", -1},
    {"feature_truncation", -1},
    {"intergenic_variant", -1},
};


// constructor
BlockSummary::BlockSummary(const char *chrnm, int num, int start, int size,
             int chrbeg, int chrend, const char *pat):
          blockIdx(num), blockStart(start), blockSize(size),
          chrBegin(chrbeg), chrEnd(chrend), isIgnored(false),
          FStat(INFINITY), pvalue(1.0), FDR(1.0), effect(0.0),
          relFStat(INFINITY), relPvalue(1.0), relFDR(1.0), relIgnore(true),
          relReject(false), numHaplo(-1), numInteresting(-1)
{
    this->chrName = strdup(chrnm);
    this->pattern = strdup(pat);
    this->numStrains = strlen(pat);
    this->numHaplo = this->numHaplotypes();
    this->makePatternUnprintable(); // this is usefull for downstream ANOVA analysis, once called, strlen(pat) will fail

}

BlockSummary::~BlockSummary()
{
    /// FIXME: need to free memory if called strdup()
    free(this->pattern);
    free(this->chrName);
}
int BlockSummary::numHaplotypes()
    {
        int numHap = -1;
        for (int str1 = 0; str1 < this->numStrains; str1++)
        {
            int hap = this->pattern[str1];
            if (this->pattern[str1] != '?' && numHap < hap)
            {
                numHap = hap;
            }
        }
        return numHap + 1;
    }

// convert digits 0-9 in string to \000..\011 (but leave '?' printable).
/// Mark: make ascii starts from 0, not '0'.
void BlockSummary::makePatternUnprintable()
{
    char *p = this->pattern;
    while (*p != 0)
    {
        if (*p != '?')
        {
            *p -= '0';
        }
        p++;
    }
}
void BlockSummary::showPatten()
{
  for (int i = 0; i < this->numStrains; i++)
  {
    char c = pattern[i];
    if ('?' == c)
    {
      std::cout << c;
    }
    else
    {
      std::cout << (char)(c + '0');
    }
  }
}

char * BlockSummary::getPatternPrintable()
{
    /// Memory Leak, be careful
    char *newpat = (char *)malloc(numStrains);
    std::memcpy(newpat, this->pattern, numStrains);

    char *p = newpat;
    while (*p != 0)
    {
        if (*p != '?')
        {
            *p += '0';
        }
        p++;
    }
    return newpat;
}

void BlockSummary::showIsCoding()  
{
    for (std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
    {
      std::cout << "\t" << (*giit).first << "\t" << (*giit).second;
    }
    std::cout << std::endl;
}

bool BlockSummary::isNumber(const std::string & str) {
{
    for (char const &c : str) 
    {
        if (std::isdigit(c) == 0) return false;
    }
    return true;
}
}
// summary of a block, read for the file.
int BlockSummary::updateCodonScore(std::string str)
  {
    /// if input strings from ANNOVAR or Ensembl-vep (INTRONIC, intergenic, SPLITE_SITE, 5PRIME_UTR, 3PRIME_UTR ...)
    /// NON_SYNONYMOUS_CODING(1), SYNONYMOUS_CODING (0), STOP have been reformat to <->
    /// Make the codon flag compatible with the old, annovar and vep annotation system
    int count = 0;
    int synonymous_count = 0;
    size_t pos = 0;
    size_t endpos = 0;
    /// missense, synonymous, stop
    /// iter through a string like: TCT/S<->CCT/P!TCT/S<->TCT/S!GAC/D<->GAC/D!GAC/D<->GAT/D
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
        int X = 'X' - 'A'; // 'X' termination condon
        if (aa1 == X || aa2 == X)
        { // stop codon
          count += 2;
          this->numInteresting = count;
          return 3; //"stop_codon";
        }
      }
      else
      {
          synonymous_count++; 
      }
      /// FIXED: iterative find < >
      pos = endpos; // find 
    }
    /// splice_donor, or splice_acceptor
    this->numInteresting = count;
    if (str.find("SPLICE_SITE") != std::string::npos)
    {
      return 2; //splicing";
    }
    
    if (count > 0)
    {
      return 1; //codon_change";
    }

    if ((count == 0) && (synonymous_count > 0))
    {
        return 0; // synonymous
    }

    /// if input string from Ensemble VEP Impact (HIGH, MODERATE, LOW, MODIFIER) or consequence string
    /// see details here: https://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html

    if (str.find("HIGH") != std::string::npos)
    {
        return 2; // Frameshift, stop gain, stop lost ...
    }
    else if (str.find("MODERATE") != std::string::npos)
    {
        return 1; // 	Missense
    }
    else if (str.find("LOW") != std::string::npos)
    {
        return 0; // Synonymous
    }
    else if (str.find("MODIFIER") != std::string::npos)
    {
        return -1; // 5 prime UTR, TF binding site varian, Intergenic ...
    } 
    else if (CSQs.find(str) != CSQs.end())
    {
        // if str is a consequence string, return the correspo
        return CSQs[str]; 
    }

    return -1; //no_codon_change"; include INTRONIC, intergenic, 5PRIME_UTR, 3PRIME_UTR
  }

void BlockSummary::updateGeneIsInteresting()
{ // BY
    bool keep = false;
    for (std::map<std::string, std::string>::iterator giit = this->geneIsCodingMap.begin(); giit != this->geneIsCodingMap.end(); giit++)
    {
        if (this->isNumber(giit->second)) 
        {
        // Our old SNP annotation only use NCBI database, which use condon flag [0,1,2,3] to indicate coding change
        // since we now using ANNOVAR to annotate SNPs, there's no 0,1s any more in the haploblock output files
        // we only keep this line here for compatibe issues
        this->geneIsInteresting[giit->first] =  std::stoi(giit->second);  
        } 
        else
        { 
            this->geneIsInteresting[giit->first] = this->updateCodonScore(giit->second);
        }

        if (this->geneIsInteresting[giit->first] >= 0 )
            keep |= true;
    }
    this->isIgnored = !keep;
}

// Return true if pval is above cutoff or FStat is below it.
bool BlockSummary::isCutoff(bool isCategorical, float cutoff)
{
    return (isCategorical) ? (this->FStat < cutoff) : (this->pvalue > cutoff);
}



BlocksComparator::BlocksComparator(bool isCat) : isCategorical(isCat){}

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
  }


GeneSummary::GeneSummary(bool ignoreDefault) : isIgnored(ignoreDefault){}
GeneSummary::~GeneSummary() {}


// Genes comparator -- compares by best blocks in gene.
GenesComparator::GenesComparator(bool isCat) : bcomp(BlocksComparator(isCat)){};
bool GenesComparator::operator()(const GeneSummary *pg1, const GeneSummary *pg2) const
{
    BlockSummary *pb1 = pg1->blocks[0];
    BlockSummary *pb2 = pg2->blocks[0];
    return bcomp(pb1, pb2);
}


GhmapWriter::GhmapWriter(char *outputFileName, char *datasetName, bool categorical, float pvalCutoff, const char* datasetFormat):
_dataset_name(datasetName), _dataset_format(datasetFormat), isCategorical(categorical), pvalueCutoff(pvalCutoff)
{
    os = std::ofstream(outputFileName);
    if (!os.is_open())
    {
        std::cout << "Open of file \"" << outputFileName << "\" failed: ";
        std::exit(1);
    }
    os << "##" << _dataset_name <<"\t"<<"FORMAT="<< _dataset_format << std::endl;
    // use datasetname to specify the codon flag explicity
    std::string dn(_dataset_name);
    upcase(dn);
    if  ((dn.find("INDEL") != std::string::npos) || (dn.find("_SV") != std::string::npos))
    {
        os << "##CodonFlag\t0:Low\t1:Moderate\t2:High\t-1:Modifier"<<std::endl;
    }
    else 
    {
    // Non-Coding -> (INTRONIC,intergenic,5PRIME_UTR,3PRIME_UTR)
        os << "##CodonFlag\t0:Synonymous\t1:Non-Synonymous\t2:Splicing\t3:Stop\t-1:Non-Coding"<<std::endl;
    }
}


GhmapWriter::~GhmapWriter() 
{
    this->os.close();
}

void GhmapWriter::sortStrainsByPheno(std::vector<std::vector<float>> &phenvec, std::vector<int> &strOrderVec)
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


// renumber eqclasses in pattern so the increase from left to right
void GhmapWriter::writeSortedPattern(char *pattern, std::vector<int> &strOrderVec)
{
    /// int numHaplo = numHaplotypes(pattern);
    int numStrains = strOrderVec.size();
    int numHaplo = -1;
    for (int str1 = 0; str1 < numStrains; str1++)
    {
        int hap = pattern[str1];
        if (pattern[str1] != '?' && numHaplo < hap)
        {
            numHaplo = hap;
        }
    }
    numHaplo += 1;

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
            this->os << (char)(sortedEqMap[eqclass] + '0');
        }
        else
        {
            this->os << '?';
        }
    }
    // memory leak?
    free(sortedEqMap);
}

void GhmapWriter::writeStrainNameAndValue(std::vector<std::vector<float>> &phenvec, std::vector<int> &strOrderVec)
{
    // output strain names.
    std::vector<int>::iterator stoEnd1 = strOrderVec.end();
    this->os << "##";
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        this->os << strainAbbrevs.eltOf(str1);
        if (stoIt1 + 1 < stoEnd1)
        {
            this->os << "\t";
        }
    }
    this->os << std::endl;


    // output phenotype values.
    this->os << "##";
    for (std::vector<int>::iterator stoIt1 = strOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        int str1 = *stoIt1;
        if (phenvec[str1].size() > 1) {
            for (unsigned i = 0; i < phenvec[str1].size() - 1; i++)
                this->os << phenvec[str1][i] << ",";
        }
        this->os << phenvec[str1].back();

        if (stoIt1 + 1 < stoEnd1)
        {
            this->os << "\t";
        }
    }
    this->os << std::endl;
}

void GhmapWriter::writeExpressionNames(std::vector<std::string> &exprOrderVec)
{
    this->os << "##GeneExprMapOrder\t";
    std::vector<std::string>::iterator stoEnd1 = exprOrderVec.end();
    for (auto stoIt1 = exprOrderVec.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        this->os << *stoIt1;
        if (stoIt1 + 1 < stoEnd1)
        {
            this->os << ";";
        }
    }
    this->os << std::endl;
}

void GhmapWriter::writeHeaders(std::vector<std::string> &header)
{
    this->os << "#";
    std::vector<std::string>::iterator stoEnd1 = header.end();
    for (auto stoIt1 = header.begin(); stoIt1 != stoEnd1; stoIt1++)
    {
        this->os << *stoIt1;
        if (stoIt1 + 1 < stoEnd1)
        {
            this->os << "\t";
        }
    }
    this->os << std::endl;
}

/// -k
/// print a line of the blocks file.
void GhmapWriter::showBlockSum(BlockSummary *pb, std::vector<int> &strOrderVec)
{
    os << pb->blockIdx;
    if (pb->isIgnored)
    {
        os << "\t(IGNORED)";
    }
    else 
    {
        os <<"\t";
    }
    os << "\t" << pb->blockStart << "\t" << pb->blockSize
        << "\t" << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t";

    this->writeSortedPattern(pb->pattern, strOrderVec);

    os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) <<"\t" << pb->effect << "\t" <<pb->FDR;
    if (!pb->relIgnore) os << "\t" <<pb->relPvalue<<"\t"<<pb->relFDR; 

    // gene names and coding bits
    std::map<std::string, std::string> &geneIsCodingMap = pb->geneIsCodingMap;
    for (std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
    { ///EDITED BY
        os << "\t" << (*giit).first << "\t" << pb->geneIsInteresting[(*giit).first];
    }
    os << std::endl;
}

/// -a
/// BY addition, output as such:
/// write all blocks that overlapped a gene
void GhmapWriter::showGeneBlockByBlock(BlockSummary *pb, std::vector<int> &strOrderVec)
  {
    std::map<std::string, std::string> &geneIsCodingMap = pb->geneIsCodingMap;
    std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin();
    if (giit == geneIsCodingMap.end())
    { // no genes in this block
      os << "None\tNone\t";
      this->writeSortedPattern(pb->pattern, strOrderVec);
      os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) <<"\t"<< pb->effect << "\t" << pb->FDR;
      if (!pb->relIgnore) os << "\t" << pb->relPvalue<<"\t"<<pb->relFDR; 
      os << "\t" << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t"
         << pb->blockIdx << "\t"<< pb->blockStart << "\t" << pb->blockSize << "\t-------------" << std::endl;
    }
    else
    {
      for (std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin(); giit != geneIsCodingMap.end(); giit++)
      {
        os << (*giit).first << "\t" << pb->geneIsInteresting[(*giit).first] << "\t";
        this->writeSortedPattern(pb->pattern, strOrderVec);
        os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) <<"\t"<< pb->effect  << "\t" <<pb->FDR;
        if (!pb->relIgnore)  os << "\t" << pb->relPvalue<<"\t"<<pb->relFDR; 
        os << "\t" <<pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd
           << "\t" <<pb->blockIdx << "\t"<< pb->blockStart << "\t" << pb->blockSize;
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

/// -k
void GhmapWriter::showBlockSums(std::vector<BlockSummary *> &blocks, float cutoff, std::vector<int> &strOrderVec)
{
    for (unsigned i = 0; i < blocks.size(); i++)
    {
      if (!blocks[i]->isCutoff(isCategorical, cutoff))
      {
        this->showBlockSum(blocks[i], strOrderVec);
      }
    }
}

void GhmapWriter::showGeneBlockByBlocks(std::vector<BlockSummary *> &blocks, float cutoff, std::vector<int> &strOrderVec, bool filterCoding)
{
    for (unsigned i = 0; i < blocks.size(); i++)
    {
        if (blocks[i]->isIgnored && filterCoding)
        {
            continue;
        }
        if (!blocks[i]->isCutoff(isCategorical, cutoff))
        {
            this->showGeneBlockByBlock(blocks[i], strOrderVec);
        }
    }
}

/// -m
/// write all blocks that overlapped to a gene
/// similar to -a
void GhmapWriter::showGeneBlockSums(std::vector<BlockSummary *> &blocks,
                       float cutoff, std::vector<int> &strOrderVec, std::vector<GeneSummary *> genesList)
{
    std::vector<std::string> genesOver; // all genes that have been already printed

    /// iter all blocks
    for (unsigned i = 0; i < blocks.size(); i++)
    {
        if (!blocks[i]->isCutoff(isCategorical, cutoff))
        {
            BlockSummary *pb = blocks[i];
            std::map<std::string, std::string> &geneIsCodingMap = pb->geneIsCodingMap;
            std::map<std::string, std::string>::iterator giit = geneIsCodingMap.begin();
            if (geneIsCodingMap.size() == 0)
            { // no genes in this block
                os << "None\tNone\t";
                this->writeSortedPattern(pb->pattern, strOrderVec);
                os << "\t" << (isCategorical ? pb->FStat : pb->pvalue) << "\t" << pb->effect << "\t";
                os << pb->chrName << "\t" << pb->chrBegin << "\t" << pb->chrEnd << "\t";
                os << pb->blockIdx << "\t"<< pb->blockStart<< "\t" << pb->blockSize << std::endl;
            }
            else
            {
                /// iter block's overlapped gene
                for (; giit != geneIsCodingMap.end(); giit++)
                {
                    /// if already printed, skip
                    std::string gname = giit->first;
                    if (std::find(genesOver.begin(), genesOver.end(), gname) != genesOver.end())
                    {
                        continue; 
                    }
                    genesOver.push_back(gname);
                    /// if not printed, find the gene in GeneSummary, write all blocks overlapped that gene
                    for (std::vector<GeneSummary *>::iterator git = genesList.begin(); git != genesList.end(); git++)
                    {
                        if ((*git)->name == gname) // only write selected gene
                        {
                            // write all blocks taht overlapped a gene
                            for (unsigned j = 0; j < ((*git)->blocks).size(); j++)
                            {
                                BlockSummary *pb_t = (*git)->blocks[j];
                                os << gname << "\t" << pb_t->geneIsInteresting[(*giit).first] << "\t";
                                this->writeSortedPattern(pb_t->pattern, strOrderVec);
                                os << "\t" << (isCategorical ? pb_t->FStat : pb_t->pvalue) <<"\t"<< pb_t->effect;
                                os << "\t" << pb_t->FDR;
                                if (!pb_t->relIgnore) os << "\t" << pb_t->relPvalue<<"\t"<<pb_t->relFDR <<"\t"; 
                                os << pb_t->chrName << "\t" << pb_t->chrBegin << "\t" << pb_t->chrEnd << "\t";
                                os << pb_t->blockIdx << "\t"<< pb_t->blockStart<< "\t" << pb_t->blockSize;
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

/// aggregate all overlapped blocks by finding the best pvalue block.
/// codon_flag indicator is >= 1 if any interesting coding change exist
/// only write the coordinate of the best block.
void GhmapWriter::showGeneBestBlockSums(std::vector<GeneSummary *> geneList, float cutoff,
                                        std::vector<int> &strOrderVec,  bool filterCoding)
{
    for (std::vector<GeneSummary *>::iterator git = geneList.begin(); git != geneList.end(); git++)
    {

        if ((*git)->isIgnored)
            continue; //simply don't write these genes
        BlockSummary *pBestBlock = (*git)->blocks[0]; // only write the first best block after GenesComparator
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
        bool isCoding = false;
        bool hasInteresting = false;
        bool hasSpliceChange = false;
        //ofstream debug_log;
        //debug_log.open("debug.log",ios::app);
        // iter all blocks that overlap a gene, find the interesting change, summaried gene blocks
        for (std::vector<BlockSummary *>::iterator blit = (*git)->blocks.begin(); blit != (*git)->blocks.end(); blit++)
        {
            if ((*blit)->isCutoff(isCategorical, cutoff))
            {
                break;
            }
            else
            {
                //if the gene had SNPs marked as NON_SYNONYMOUS_CODING, with <->, or as SPLICE_SITE, isCoding is true
                isCoding |= ((*blit)->geneIsInteresting[gname] >=0 ); // synousmous or missene
                hasInteresting |= ((*blit)->numInteresting > 0);   //has a major amino acid change, missense 
                // hasSpliceChange |= (((*blit)->geneIsCodingMap[gname]).find("SPLICE_SITE") != std::string::npos);
                // need to compatible with High impact or SPLITE_SITE 
                hasSpliceChange |= ((*blit)->geneIsInteresting[gname] >=2 );
            }
        }

        //New thing: -1 means not coding, 0 means coding but not important
        //Anything else is the number of important
        int codingCode = -1;
        if(isCoding) codingCode = 0;
        if(hasInteresting) codingCode = 1;
        if (hasSpliceChange) codingCode = 2;

        if ( (isCoding || !filterCoding) || hasSpliceChange) 
        {
            
            // os << gname << "\t" << pBestBlock->geneIsInteresting[(*git)->name] << "\t";
            /// for debug
            // os << gname << "\t" << pBestBlock->geneIsInteresting[(*git)->name] << "\t" << pBestBlock->geneIsCodingMap[gname] << "\t";
            os << gname << "\t" << codingCode << "\t";
            this->writeSortedPattern(pBestBlock->pattern, strOrderVec);
            os << "\t" << (isCategorical ? pBestBlock->FStat : pBestBlock->pvalue);
            os << "\t" << pBestBlock->effect <<"\t"<< pBestBlock->FDR;
            if (!pBestBlock->relIgnore) os << "\t" << pBestBlock->relPvalue<<"\t"<< pBestBlock->relFDR; 
            os << "\t" << pBestBlock->chrName << "\t" << pBestBlock->chrBegin << "\t" << pBestBlock->chrEnd;
            //os << "\t" << pBestBlock->blockIdx << "\t" << pBestBlock->blockStart << "\t" << pBestBlock->blockSize;

            // write gene expression values
            if (geneExprMap.find(ugname) == geneExprMap.end())
            {
                // no data for that gene name
                os << "\t-------------" << std::endl;
            }
            else
            {
                os << "\t" << geneExprMap[ugname] << std::endl;
            }
        }
    }
}


void readBlockSummary(char *fname)
{
    ColumnReader rdr(fname, (char *)"\t");

    int numtoks;
    while ((numtoks = rdr.getLine()) >= 0)
    {
        if (rdr.getCurrentLineNum() < 1) continue;
        // file has "Chromosome\tBlockNum\tBlockstartSNPnum\tSize\tChrbegin\tChrend\tPattern\tgene1\tcodon1\t...\n"
        BlockSummary *pBlock = new BlockSummary(rdr.getToken(0).c_str(),
                                                std::stoi(rdr.getToken(1)),
                                                std::stoi(rdr.getToken(2)),
                                                std::stoi(rdr.getToken(3)),
                                                std::stoi(rdr.getToken(4)),
                                                std::stoi(rdr.getToken(5)),
                                                rdr.getToken(6).c_str());

        // For each gene, add gene summary to the gene table (if not there already),
        // and insert block in the block set.
        // This iterates over list of alternating gene name/codon tokens.
        for (std::vector<std::string>::iterator git = rdr.begin() + 7; git != rdr.end(); git += 2)
        {
            // *(git+1) is codon for gene name (*git).  Convert from string to bool using cond expr.
            pBlock->geneIsCodingMap[*git] = *(git + 1); // e.g. Oct4: TGA/X<->GGA/R
            std::unordered_map<std::string, GeneSummary *>::iterator entIt = geneTable.find(*git);

            if (entIt == geneTable.end())
            {
                // create new entry.
                // cout << "  Creating new gene summary for " << *git << endl;
                geneTable[*git] = new GeneSummary(false);
                geneTable[*git]->name = *git;
            }
            // Add block to GeneSummary
            //      cout << "  Adding block for " << *git << ": " << pLastBlock->blockIdx << endl;
            geneTable[*git]->blocks.push_back(pBlock);
        }
        pBlock->updateGeneIsInteresting(); //
        // Add to blocks.
        blocks.push_back(pBlock);
        // cout << "Gene names for block " << blocks.back().blockIdx << ": " << blocks.back().geneNames << endl;
        
    }
    // Separate pass to fix up patterns.
    // patterns are all the same length.
    // FIXME: this causes a seg fault when there are no blocks (happens under funny conditions).
    if (blocks.empty())
    {
        std::cerr<<"Input filename has not content : "<<fname<<std::endl;
    }
    numStrains = blocks[0]->numStrains; // this line is required
}


/// -k 
void writeBlockSums(bool isCategorical, char *outputFileName,
                    char *datasetName, std::vector<std::vector<float> > &phenvec,
                    std::vector<BlockSummary *> &blocks, float pvalueCutoff)
{
    int numStrains = phenvec.size();
    std::vector<int> strOrderVec(numStrains);
    GhmapWriter writer(outputFileName, datasetName, isCategorical, pvalueCutoff, "block-oriented");
    writer.sortStrainsByPheno(phenvec, strOrderVec);
    writer.writeExpressionNames(geneExprHeader);
    writer.writeStrainNameAndValue(phenvec, strOrderVec);
    writer.os << "#BlockIdx\tIGNORED\tBlockStart\tBlockSize\tChr\tChrStart\tChrEnd\tHaplotype\t";
    writer.os << (isCategorical ? "FStat" : "Pvalue");
    writer.os << "\tEffectSize\tFDR\t";
    if (!blocks[0]->relIgnore) writer.os << "PopPvalue\tPopFDR\t";
    writer.os <<"GeneName\tCodonFlag\n";
    writer.showBlockSums(blocks, pvalueCutoff, strOrderVec);
}

/// -a
void writeGeneBlockByBlocks(bool isCategorical, char *outputFileName, char *datasetName,
                            std::vector<std::vector<float>> &phenvec,std::vector<BlockSummary *> &blocks, float cutoff, bool filterCoding)
{
    int numStrains = phenvec.size();
    std::vector<int> strOrderVec(numStrains); // will contain strain indices.
    GhmapWriter writer(outputFileName, datasetName, isCategorical, cutoff, "gene-oriented");
    writer.sortStrainsByPheno(phenvec, strOrderVec);
    writer.writeExpressionNames(geneExprHeader);
    writer.writeStrainNameAndValue(phenvec, strOrderVec);
    //writer.writeHeaders(header);
    writer.os << "#GeneName\tCodonFlag\tHaplotype\t";
    writer.os << (isCategorical ? "FStat" : "Pvalue");
    writer.os << "\tEffectSize\tFDR\t";
    if (!blocks[0]->relIgnore) writer.os << "PopPvalue\tPopFDR\t";
    writer.os << "Chr\tChrStart\tChrEnd\t";
    writer.os << "BlockIdx\tBlockStart\tBlockSize\tGeneExprMap\n";
    writer.showGeneBlockByBlocks(blocks, cutoff, strOrderVec, filterCoding);
}

/// -m
void writeGeneBlockSums(bool isCategorical, char *outputFileName, char *datasetName,
                        std::vector<std::vector<float>> &phenvec, std::vector<BlockSummary *> &blocks, float pvalueCutoff)
{
    int numStrains = phenvec.size();
    std::vector<int> strOrderVec(numStrains); // will contain strain indices.
    GhmapWriter writer(outputFileName, datasetName, isCategorical, pvalueCutoff,"gene-oriented");
    writer.sortStrainsByPheno(phenvec, strOrderVec);
    writer.writeExpressionNames(geneExprHeader);
    writer.writeStrainNameAndValue(phenvec, strOrderVec);
    writer.os << "#GeneName\tCodonFlag\tHaplotype\t";
    writer.os << (isCategorical ? "FStat" : "Pvalue");
    writer.os << "\tEffectSize\tFDR\t";
    if (!blocks[0]->relIgnore) writer.os << "PopPvalue\tPopFDR\t";
    writer.os << "Chr\tChrStart\tChrEnd\t";
    writer.os << "BlockIdx\tBlockStart\tBlockSize\tGeneExprMap\n";

    GenesComparator gcomp(isCategorical);
    std::vector<GeneSummary *> genes;
    genes.reserve(geneTable.size());
    //transform(geneTable.begin(), geneTable.end(), back_inserter(genes), select2nd<std::unordered_map<string, GeneSummary *> >());
    std::transform(geneTable.begin(), geneTable.end(), std::back_inserter(genes),
              std::bind(&std::unordered_map<std::string, GeneSummary *>::value_type::second, std::placeholders::_1 ));

    std::sort(genes.begin(), genes.end(), gcomp);
    writer.showGeneBlockSums(blocks, pvalueCutoff, strOrderVec, genes);
}

// Write gene-oriented summary. default
void writeGeneSums(bool isCategorical, char *outputFileName,
                   char *datasetName, std::vector<std::vector<float> > &phenvec,
                   std::vector<BlockSummary *> &blocks, float pvalueCutoff, bool filterCoding)
{
    int numStrains = phenvec.size();
    std::vector<int> strOrderVec(numStrains);
    const char * fmt = "gene-summaried; To get exact CodonFlag coordinates, re-run ghmap with -k/m/a";
    GhmapWriter writer(outputFileName, datasetName, isCategorical, pvalueCutoff, fmt);
    writer.sortStrainsByPheno(phenvec, strOrderVec);
    writer.writeExpressionNames(geneExprHeader);
    writer.writeStrainNameAndValue(phenvec, strOrderVec);
    // header
    writer.os << "#GeneName\tCodonFlag\tHaplotype\t";
    writer.os << (isCategorical ? "FStat" : "Pvalue");
    writer.os << "\tEffectSize\tFDR\t";
    if (!blocks[0]->relIgnore) writer.os << "PopPvalue\tPopFDR\t";
    writer.os << "Chr\tChrStart\tChrEnd\tGeneExprMap\n";

    GenesComparator gcomp(isCategorical);
    // Copy genesTable values into a vector and sort using GenesComparator
    std::vector<GeneSummary *> geneList;
    geneList.reserve(geneTable.size());
    std::transform(geneTable.begin(), geneTable.end(), std::back_inserter(geneList),
              std::bind(&std::unordered_map<std::string, GeneSummary *>::value_type::second, std::placeholders::_1 ));
    std::sort(geneList.begin(), geneList.end(), gcomp); // Compare blocks
    writer.showGeneBestBlockSums(geneList, pvalueCutoff, strOrderVec,  filterCoding);
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
        if (rdr.getCurrentLineNum() < 1) continue;
        // A typical line: "Myc	PAPAAAAAAAPAA"
        std::string geneName = rdr.getToken(0);
        std::string present = rdr.getToken(1);

        upcase(geneName);

        // Gene names in expr file are all upper case, so use upper case comparisons.
        // Don't really need to do upcase here, but doing it anyway in case someone
        // substitutes a file with mixed-case names.
        geneExprMap[geneName] = present;
    }
    geneExprHeader = rdr.getHeaderLines().back(); // only get last header vector
}

// Read a file of quantitative phenotypes.
void readQPhenotypes(char *fname, std::vector<std::vector<float>> &phenvec, Dynum<std::string> & strainAbbrevs)
{
    ColumnReader rdr(fname, (char *)"\t");

    int numtoks;
    int qtok = 1; // floating point value token
    std::string::size_type sz; 
    while ((numtoks = rdr.getLine()) >= 0)
    {
        if (rdr.getCurrentLineNum() < 1) continue;
        // file has "Abbrev\tFullname\tValue1,Value2,Value3\n"

        std::string strain_abbrev = rdr.getToken(0);
        //qphen.push_back(std::stof(rdr.getToken(1)));
        //    int strIdx = strainAbbrevs.hasIndex(strain_abbrev);
        int strIdx = strainAbbrevs.addElementIfNew(strain_abbrev);
        if (strIdx < 0)
        {
           std::cout << "Undefined strain abbrev: " << strain_abbrev << std::endl;
        }

        // test whether it's 2 columns or 3 columns format
        try 
        { 
            std::stof(rdr.getToken(qtok), &sz); 
            //qtok = 1;
        } 
        catch(std::exception& ia) 
        { 
            qtok ++;
        } 

        // now parse the floating point values
        std::vector<std::string> qphen;
        rdr.split(rdr.getToken(qtok), (char *)",", qphen);
        for (auto &q: qphen)
        {
            if (q.empty()) continue;
            phenvec[strIdx].push_back(std::stof(q));
        }
    }
}

// Read a file of categorical phenotypes.
void readCPhenotypes(char *fname, std::vector<std::vector<float>> &phenvec, Dynum<std::string> & strainAbbrevs)
{

    ColumnReader rdr(fname, (char *)"\t");
    Dynum<std::string> categories; // assigned the distinct categories consecutive indices, starting at 0.

    catNames.resize(numStrains);

    int numtoks;
    int qtok = 1; // floating point value token
    std::string::size_type sz; 
    while ((numtoks = rdr.getLine()) >= 0)
    {
        if (rdr.getCurrentLineNum() < 1) continue;
        // file has "Abbrev\t\Value"
        if (numtoks != 2)
        {
            std::cout << "Warning: numtoks = " << numtoks << std::endl;
        }
        // test whether it's 2 columns or 3 columns format
        try 
        { 
            std::stof(rdr.getToken(qtok), &sz); 
            //qtok = 1;
        } 
        catch(std::exception& ia) 
        { 
            qtok ++;
        } 

        // FIXME: some unnecessary string copies
        std::string strain_abbrev = rdr.getToken(0);
        std::string catname = rdr.getToken(qtok);
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
        if (rdr.getCurrentLineNum() < 1) continue;
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
        if (rdr.getCurrentLineNum() < 1) continue;
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
        if (rdr.getCurrentLineNum() < 1) continue;
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
