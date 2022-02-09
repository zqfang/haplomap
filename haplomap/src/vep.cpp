//
// Created by Zhuoqing Fang on 2/1/22.
//

#include "vep.h"

Key::Key(std::string key) 
{
    size_t tokstart = 0;
    size_t dpos = 0;
    size_t sz = key.size();

    // split consequence strings
    std::vector<std::string> tmp;
    tokstart = 0;
    while (dpos < key.size())
    {
        dpos = key.find_first_of("_", tokstart);
        // this crock necessitated by weird segfault.
        if (dpos > sz)
        {
        dpos = sz;
        }
        std::string tok = key.substr(tokstart, dpos - tokstart);
        if (tok.size() > 0) tmp.push_back(tok);
        tokstart = dpos + 1; // don't care if it overflows
    }
    this->chrom = tmp[1];
    this->start = std::stoi(tmp[2]);
    if (tmp.size() >= 4) 
        this->end = std::stoi(tmp[3]);
    else
        this->end = this->start;
}

// constructor
VEPSummary::VEPSummary(
    std::string uploaded_variant, std::string loc, std::string seq, std::string gene,
    std::string transcript, std::string feature_type, std::string csq, std::string aa_pos, std::string aa,
    std::string codon, std::string existing_variation, std::string ind, std::string zyg, std::string imp,
    std::string var_class, std::string gene_name, std::string bt, std::string hgvsp, std::string hgvsc) : 
    idx(uploaded_variant), location(loc), allel(seq), geneid(gene), transxid(transcript),
    transtype(feature_type), protein_position(aa_pos),
    zygosity(zyg), variant_class(var_class), symbol(gene_name), biotype(bt)
{
    size_t tokstart = 0;
    size_t dpos = 0;
    size_t sz = location.size();
    size_t tok1 = location.find_first_of(":", tokstart);
    size_t tok2 = location.find_first_of("-", tokstart);

    this->chrom = location.substr(tokstart, tok1);

    /// For ensembl-vep results, the coordinates (chrStart) of
    /// indel, deletion need to -1 to get original position in vcf.
    /// SNV, insertion stay the same to position in vcf
    /// The quick trick to handle these cases is whether the location string contains "-"
    if (tok2 == std::string::npos) // no match found
    {
        this->start = std::stoi(location.substr(tok1 + 1, sz-tok1)); // (pos, len)
        this->end = this->start; 
        if (var_class == "indel") this->start --; // indel (dup) with "chr:start" format need to minus 1 
    }
    else
    {
        // NOTE: need to minius 1, since VEP made pos+1 in their annotatoin
        this->start = std::stoi(location.substr(tok1 + 1, tok2 - tok1)) - 1; 
        if (var_class == "insertion") this->start ++; // only insertion case are stay the same pos as original vcf
        this->end = std::stoi(location.substr(tok2 + 1, sz - tok2)); // empty string if snp
    }
    // split consequence strings
    tokstart = 0;
    while (dpos < csq.size())
    {
        dpos = csq.find_first_of(",", tokstart);
        // this crock necessitated by weird segfault.
        size_t csz = csq.size();
        if (dpos > csz)
        {
        dpos = csz;
        }
        std::string tok = csq.substr(tokstart, dpos - tokstart);
        consequence.addElementIfNew(tok);
        tokstart = dpos + 1; // don't care if it overflows
    }
    amino_acids.addElementIfNew(aa);
    codons.addElementIfNew(codon);
    dbsnpid.addElementIfNew(existing_variation);
    samples.addElementIfNew(ind);
    impact.addElementIfNew(imp);
    HGVSc.addElementIfNew(hgvsc);
    HGVSp.addElementIfNew(hgvsp);
}

VarirantEeffectPredictor::VarirantEeffectPredictor(char *inVEPName, char *inStrainName) : numtoks(-1), hasIND(true)
{
    CODONs = {{"TTT", "F"}, {"TTC", "F"}, {"TCT", "S"}, {"TCC", "S"}, {"TAT", "Y"}, {"TAC", "Y"}, {"TGT", "C"}, {"TGC", "C"}, 
          {"TTA", "L"}, {"TCA", "S"}, {"TAA", "X"}, {"TGA", "X"}, {"TTG", "L"}, {"TCG", "S"}, {"TAG", "X"}, {"TGG", "W"}, 
          {"CTT", "L"}, {"CTC", "L"}, {"CCT", "P"}, {"CCC", "P"}, {"CAT", "H"}, {"CAC", "H"}, {"CGT", "R"}, {"CGC", "R"}, 
          {"CTA", "L"}, {"CTG", "L"}, {"CCA", "P"}, {"CCG", "P"}, {"CAA", "Q"}, {"CAG", "Q"}, {"CGA", "R"}, {"CGG", "R"}, {"ATT", "I"}, 
          {"ATC", "I"}, {"ACT", "T"}, {"ACC", "T"}, {"AAT", "N"}, {"AAC", "N"}, {"AGT", "S"}, {"AGC", "S"}, {"ATA", "I"}, {"ACA", "T"}, 
          {"AAA", "K"}, {"AGA", "R"}, {"ATG", "M"}, {"ACG", "T"}, {"AAG", "K"}, {"AGG", "R"}, {"GTT", "V"}, {"GTC", "V"}, {"GCT", "A"}, 
          {"GCC", "A"}, {"GAT", "D"}, {"GAC", "D"}, {"GGT", "G"}, {"GGC", "G"}, {"GTA", "V"}, {"GTG", "V"}, {"GCA", "A"}, {"GCG", "A"}, 
          {"GAA", "E"}, {"GAG", "E"}, {"GGA", "G"}, {"GGG", "G"}};
    PRIOR = {{"HIGH", 2}, {"MODERATE", 1}, {"LOW", 0}, {"MODIFIER", -1}};
    CSQs = { {"transcript_ablation", "HIGH"},
             {"splice_acceptor_variant", "HIGH"},
             {"splice_donor_variant", "HIGH"},
             {"stop_gained", "HIGH"},
             {"frameshift_variant", "HIGH"},
             {"stop_lost", "HIGH"},
             {"start_lost", "HIGH"},
             {"transcript_amplification", "HIGH"},
             {"inframe_insertion", "MODERATE"},
             {"inframe_deletion", "MODERATE"},
             {"missense_variant", "MODERATE"},
             {"protein_altering_variant", "MODERATE"},
             {"splice_region_variant", "LOW"},
             {"incomplete_terminal_codon_variant", "LOW"},
             {"start_retained_variant", "LOW"},
             {"stop_retained_variant", "LOW"},
             {"synonymous_variant", "LOW"},
             {"coding_sequence_variant", "MODIFIER"},
             {"mature_miRNA_variant", "MODIFIER"},
             {"5_prime_UTR_variant", "MODIFIER"},
             {"3_prime_UTR_variant", "MODIFIER"},
             {"non_coding_transcript_exon_variant", "MODIFIER"},
             {"intron_variant", "MODIFIER"},
             {"NMD_transcript_variant", "MODIFIER"},
             {"non_coding_transcript_variant", "MODIFIER"},
             {"upstream_gene_variant", "MODIFIER"},
             {"downstream_gene_variant", "MODIFIER"},
             {"TFBS_ablation", "MODIFIER"},
             {"TFBS_amplification", "MODIFIER"},
             {"TF_binding_site_variant", "MODIFIER"},
             {"regulatory_region_ablation", "MODERATE"},
             {"regulatory_region_amplification", "MODIFIER"},
             {"feature_elongation", "MODIFIER"},
             {"regulatory_region_variant", "MODIFIER"},
             {"feature_truncation", "MODIFIER"},
             {"intergenic_variant", "MODIFIER"},
    };
    // read strain file
    if (inStrainName != nullptr)
    {
        ColumnReader sdr(inStrainName, (char *)"\t");
        while ((numtoks = sdr.getLine()) >= 0)
        {
            if (sdr.getCurrentLineNum() < 1)
                continue;
            // file has "Abbrev\tFullname\tValue1,Value2,Value3\n"
            std::string abbrev = sdr.getToken(0);
            strainAbbrevs.addElementIfNew(abbrev);
        }
    }
    // read VEP header
    this->readHeader(inVEPName, (char *)"\t");
}

VarirantEeffectPredictor::~VarirantEeffectPredictor()
{
    for (auto *p : data)
        delete p;
    
}

void::VarirantEeffectPredictor::upcase(std::string & str)
{
  // upase strings
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

void::VarirantEeffectPredictor::lowercase(std::string & str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

void VarirantEeffectPredictor::readHeader(char *inFileName, char *delemiter)
{
    ColumnReader rdr(inFileName, delemiter);
    // read header first to get column names
    while ((numtoks = rdr.getLine()) >= 0)
    {
        if (rdr.getCurrentLineNum() >= 1)
            break; // read headers
    }
    std::vector<std::string> _columns = rdr.getHeaderLines().back();
    for (unsigned i = 0; i < _columns.size(); i ++)
    {
        columns[_columns[i]] = i;
    }
    // wheter has IND column
    if (columns.find("IND") == columns.end())
    {
        this->hasIND = false;
    }

}
/// varType (variant classes): https://m.ensembl.org/info/genome/variation/prediction/classification.html
void VarirantEeffectPredictor::readVEP(char *inVEPName, char *delemiter, char* varType)
{

    ColumnReader rdr(inVEPName, delemiter);
    // read header first to get column names
    std::string location = "";
    std::string key = "SNP_";
    std::string transcript_id = "ENT";
    std::string _varType(varType);

    lowercase(_varType);
    while ((numtoks = rdr.getLine()) >= 0)
    {
        // skip header
        if (rdr.getCurrentLineNum() < 1)
            continue;
        // 
        std::string _varclass = rdr.getToken(columns["VARIANT_CLASS"]);
        lowercase(_varclass);
        // if varType is "sv", "all", write all
        if ((_varType == "snv" || _varType == "indel") && (_varclass.find(_varType) == std::string::npos)) //
            continue;
        /// if the record is no in the queried input strains, skip.
        /// this ensure that we could identify impactful variants only from the queried strains.
        if ((strainAbbrevs.size() > 0) && hasIND) // if no IND column, then we need to parse the record
        {
            if (strainAbbrevs.hasIndex(rdr.getToken(columns["IND"])) < 0) continue;
        }
        /// if a heterozyote is annotated, treat them as the same to the reference genome, skip here
        if (rdr.getToken(columns["ZYG"]) == "HET")
            continue;
        
        /// aggregate results groupby location and transcript
        transcript_id = rdr.getToken(columns["Feature"]);
        if ((geneCodingMap.find(key) != geneCodingMap.end()) && (geneCodingMap[key].find(transcript_id) != geneCodingMap[key].end()))
        {
            /// now add annotations
            VEPSummary *pVEP = geneCodingMap[key][transcript_id];

            // aggreate individual sample results onto transcript level
            if (hasIND) pVEP->samples.addElementIfNew(rdr.getToken(columns["IND"]));
            //pVEP->zygosity.append(","+rdr.getToken(14)); // only "HOM"
            pVEP->impact.addElementIfNew(rdr.getToken(columns["IMPACT"]));
            // if (pVEP->HGVS != rdr.getToken(32))
            //     pVEP->impact.append(","+rdr.getToken(32));
            pVEP->dbsnpid.addElementIfNew(rdr.getToken(columns["Existing_variation"]));
            pVEP->codons.addElementIfNew(rdr.getToken(columns["Codons"]));
            pVEP->amino_acids.addElementIfNew(rdr.getToken(columns["Amino_acids"]));
            //pVEP->variant_class.addElementIfNew(rdr.getToken(19));
            size_t tokstart = 0;
            size_t dpos = 0;
            std::string csq = rdr.getToken(columns["Consequence"]);
            while (dpos < csq.size())
            {
                dpos = csq.find_first_of(",", tokstart);
                // this crock necessitated by weird segfault.
                size_t csz = csq.size();
                if (dpos > csz) dpos = csz;
                pVEP->consequence.addElementIfNew(csq.substr(tokstart, dpos - tokstart));
                tokstart = dpos + 1; // don't care if it overflows
            }
        }
        else
        {
            location = rdr.getToken(columns["Location"]);
            VEPSummary *pRecord = new VEPSummary(rdr.getToken(columns["Uploaded_variation"]), // alread trim # when read header
                                                 rdr.getToken(columns["Location"]),
                                                 rdr.getToken(columns["Allele"]),
                                                 rdr.getToken(columns["Gene"]),
                                                 rdr.getToken(columns["Feature"]),
                                                 rdr.getToken(columns["Feature_type"]),
                                                 rdr.getToken(columns["Consequence"]),
                                                 rdr.getToken(columns["Protein_position"]),
                                                 rdr.getToken(columns["Amino_acids"]),
                                                 rdr.getToken(columns["Codons"]),
                                                 rdr.getToken(columns["Existing_variation"]),
                                                 rdr.getToken(columns["IND"]),
                                                 rdr.getToken(columns["ZYG"]),
                                                 rdr.getToken(columns["IMPACT"]),
                                                 rdr.getToken(columns["VARIANT_CLASS"]),
                                                 rdr.getToken(columns["SYMBOL"]),
                                                 rdr.getToken(columns["BIOTYPE"]),
                                                 rdr.getToken(columns["HGVSc"]),
                                                 rdr.getToken(columns["HGVSp"]));
            this->data.push_back(pRecord); // for deconstuctor
            /// FIXME: indel annotation vep input also have entries indel, insertion, deletion. How to handel these?
            int var_len = pRecord->end - pRecord->start;
            if (pRecord->variant_class == "SNV")
            {
                key = "SNP_" + pRecord->chrom + "_" + std::to_string(pRecord->start);
            }
            else if ((pRecord->variant_class == "indel") ||
                     ((pRecord->variant_class == "insertion") && var_len < 50) ||  
                     ((pRecord->variant_class ==  "deletion") && var_len < 50) )
            {
                /// NOTE: VEP indels are > 2bp, else it will annotate as deletions and insertions. 
                /// So, defined 1 bp del or ins as Indels for downstream analysis
                /// see docs: https://m.ensembl.org/info/genome/variation/prediction/classification.html  
                /// we force var_len < 50 bp to be indels
                key = "INDEL_" + pRecord->chrom + "_" + std::to_string(pRecord->start);
            }
            else
            {
                std::string _svtype = pRecord->variant_class.substr(0,3);
                upcase(_svtype);
                key = "SV_" + pRecord->chrom + "_" 
                            + std::to_string(pRecord->start) + "_" 
                            + std::to_string(pRecord->end) + "_" 
                            + _svtype;
            }
            this->geneCodingMap[key][transcript_id] = pRecord;
        }
    }
    std::unordered_map<std::string, std::unordered_map<std::string, VEPSummary *>>::iterator giit = geneCodingMap.begin();
    for (; giit != geneCodingMap.end(); giit++)
        this->keys.push_back(giit->first);
    std::sort(this->keys.begin(), this->keys.end(), [&](std::string a, std::string b ){return compareKey(Key(a), Key(b));});
}
std::string VarirantEeffectPredictor::codonChange(VEPSummary * pRecord)
{   // aggregate all condon changes groupby (location, transcript)
    std::string _expr;
    for (int j = 0; j < pRecord->codons.size(); j++)
    {
        std::string codon = pRecord->codons.eltOf(j);
        if (codon == "-") continue;
        // std::string _codon(codon);
        upcase(codon);
        int pos = codon.find_first_of("/", 0);
        std::string ref = codon.substr(0, pos); // keep original case sensitive strings for expr
        std::string alt = codon.substr(pos+1, codon.size() - pos);
        //std::string _ref = _codon.substr(0, pos);
        //std::string _alt = _codon.substr(pos+1, codon.size() - pos);
        std::string expr = ref+"/"+CODONs[ref]+"<->"+alt+"/"+CODONs[alt];
        if (_expr.size() > 0 )
            _expr.append("!"+expr);
        else
            _expr = expr; 
    }
    return _expr;

}
bool VarirantEeffectPredictor::compareKey(Key key1, Key key2)
{   // sort first by chrom, next by start, then by end 
    if (key1.chrom != key2.chrom)
    {
        // deal with complexities of comparing chr names.
        // put "M" last
        if (key2.chrom == "M")
        {
        return true;
        }
        else if (key1.chrom == "M")
        {
        return false;
        }
        // Y next to last
        else if (key2.chrom == "Y")
        {
        return true;
        }
        else if (key1.chrom == "Y")
        {
        return false;
        }
        // X just before Y
        else if (key2.chrom == "X")
        {
        return true;
        }
        else if (key1.chrom == "X")
        {
        return false;
        }
        // otherwise, in numerical order
        int numChr1 = std::stoi(key1.chrom); // char* to integer
        int numChr2 = std::stoi(key2.chrom);
        if (numChr1 < numChr2) return true;
        else if (numChr1 > numChr2) return false;
        else return false;
    }
    // chromosome names were same, so use position
    // key1=key2 for primary condition, go to secondary
    if (key1.start < key2.start) return true;
    if (key2.start < key1.start) return false;
    // key1=key2 for primary condition, go to third
    if (key1.end < key2.end) return true;
    if (key2.end < key1.end) return false;
    return false;
} 


void VarirantEeffectPredictor::writeVEPImpact(char* outFileName)
{
    std::ofstream csqos(outFileName);
    for (auto & k: this->keys)
    { // iter variant
        csqos << k;
        // aggregate impact to variant level
        std::unordered_map<std::string, VEPSummary *>::iterator transxit = geneCodingMap[k].begin();
        Dynum<std::string> csq;
        for (; transxit != geneCodingMap[k].end(); transxit++)
        { // iter transcripts
            VEPSummary *pRecord = transxit->second;
            if (pRecord->symbol == "-") continue;        
            for (int i = 0; i < pRecord->impact.size(); i++)
            {
                std::string s = pRecord->symbol +"\t" + pRecord->impact.eltOf(i);
                if (csq.hasIndex(s) < 0)
                { // remove duplicate entries
                    csqos<< "\t" << s <<"\t";
                    for (int j=0; j < pRecord->samples.size(); j ++)
                    {
                        csqos<<pRecord->samples.eltOf(j);
                        if ((j +1 ) < pRecord->samples.size()) 
                            csqos <<",";
                    }
                    csq.addElementIfNew(s);
                }
            } 
        }
        csqos << std::endl;
    }
    csqos.close();
}

void VarirantEeffectPredictor::writeVEPCsq(char* outFileName) 
{
    std::ofstream csqos(outFileName);
    for (auto & k: this->keys)
    { // iter variant
        csqos << k;
        // aggregate consequences to variant level
        std::unordered_map<std::string, VEPSummary *>::iterator transxit = geneCodingMap[k].begin();
        Dynum<std::string> csq; // aggreate transcript results to snp level
        for (; transxit != geneCodingMap[k].end(); transxit++)
        { // iter transcripts
            VEPSummary *pRecord = transxit->second;
            if (pRecord->symbol == "-") continue;        
            for (int i = 0; i < pRecord->consequence.size(); i++)
            {
                std::string s = pRecord->symbol + "\t";//+ pRecord->consequence.eltOf(i);
                std::string csq_str = pRecord->consequence.eltOf(i);
                if (CSQs[csq_str] == "HIGH")
                {
                    if (csq_str.find("splice") != std::string::npos)
                    {
                        s.append("SPLICE_SITE");
                    } else if (csq_str.find("stop") != std::string::npos)
                    {
                        s.append(this->codonChange(pRecord));
                        //s.append("STOP");
                    } else
                    {
                        s.append(csq_str);
                    }
                } 
                else if (csq_str == "missense_variant" || csq_str == "synonymous_variant")
                {
                    s.append(this->codonChange(pRecord));
                }
                else 
                {
                    //s.append(CSQs[csq_str]);
                    s.append(csq_str);
                }
                if (csq.hasIndex(s) < 0)
                { // remove duplicate entries
                    csqos<< "\t" << s ; //<<"\t"<<pRecord->samples;
                    csq.addElementIfNew(s);
                }
          
            } 
               
        }
        csqos << std::endl;
    }
    csqos.close();

}