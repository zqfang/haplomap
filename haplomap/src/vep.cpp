//
// Created by Zhuoqing Fang on 2/1/22.
//

#include "vep.h"
#include "constants.h"

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
    // Key* loc2 = Key(location);
    // this->chrom = loc2->chrom;
    // this->start = loc2->start;
    // this->end = loc2->end;
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

void VarirantEeffectPredictor::upcase(std::string & str)
{
  // upase strings
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

void VarirantEeffectPredictor::lowercase(std::string & str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

bool VarirantEeffectPredictor::isValidDNASequence(const std::string& sequence) {
    return std::all_of(sequence.begin(), sequence.end(), [](char nucleotide) {
        return nucleotide == 'a' || nucleotide == 'c' || nucleotide == 'g' || nucleotide == 't';
    });
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
    if (columns.find("VARIANT_CLASS") == columns.end())
    {
        std::cerr<<"VARIANT_CLASS column is missing, please re-run vep with `--variant_class`\n";
        std::exit(1);
    }

}

/// var_class must be lowercase
/// var_type: snv, indel, sv. 
std::string VarirantEeffectPredictor::set_key(std::string location, std::string var_class, std::string& var_type)
{
    size_t tokstart = 0;
    size_t sz = location.size();
    size_t tok1 = location.find_first_of(":", tokstart);
    size_t tok2 = location.find_first_of("-", tokstart);
    size_t tok0 = location.find_first_not_of("chr");
    int start;
    int end;
    std::string key;
    std::string chrom = location.substr(tok0, tok1);

    /// For ensembl-vep results, the coordinates (chrStart) of
    /// indel, deletion need to -1 to get original position in vcf.
    /// SNV, insertion stay the same to position in vcf
    /// The quick trick to handle these cases is whether the location string contains "-"
    if (tok2 == std::string::npos) // no match found
    {
        start = std::stoi(location.substr(tok1 + 1, sz-tok1)); // (pos, len)
        end = start; 
        if (var_class == "indel") start --; // indel (dup) with "chr:start" format need to minus 1 
    }
    else
    {
        // NOTE: need to minius 1, since VEP made pos+1 in their annotatoin
        start = std::stoi(location.substr(tok1 + 1, tok2 - tok1)) - 1; 
        if (var_class == "insertion") start ++; // only insertion case are stay the same pos as original vcf
        end = std::stoi(location.substr(tok2 + 1, sz - tok2)); // empty string if snp
    }

    int var_len = end - start;
    if (var_type == "sv") 
    {
        // for sv's inserstion, var_len is 1. 
        // set to 100, so the key is like SV_** 
        var_len = 100; 
        if (var_class == "insertion") end = start;
    }


    if (var_class == "snv")
    {
        key = "SNP_" + chrom + "_" + std::to_string(start);
    }
    else if ((var_class == "indel") ||
                ((var_class == "insertion") && var_len < 50) ||  
                ((var_class ==  "deletion") && var_len < 50) )
    {
        /// NOTE: VEP indels are > 2bp, or else it will annotate as deletions and insertions. 
        /// So, defined 1 bp del or ins as Indels for downstream analysis
        /// see docs: https://ensembl.org/info/genome/variation/prediction/classification.html
        /// we force var_len < 50 bp to be indels
        key = "INDEL_" + chrom + "_" + std::to_string(start);
    }
    else
    {
        std::string _svtype = var_class.substr(0,3);
        this->upcase(_svtype);
        key = "SV_" + chrom + "_" 
                    + std::to_string(start) + "_" 
                    + std::to_string(end) + "_" 
                    + _svtype;
    }
    return key;

}

/// varType (variant classes): 
/// https://ensembl.org/info/genome/variation/prediction/classification.html
void VarirantEeffectPredictor::readVEP(char *inVEPName, char *delemiter, char* varType)
{

    ColumnReader rdr(inVEPName, delemiter);
    // read header first to get column names
    std::string location;
    std::string key;
    std::string transcript_id = "ENT";
    std::string _varType(varType);
    std::string _varclass;

    lowercase(_varType);
    while ((numtoks = rdr.getLine()) >= 0)
    {
        // skip header
        if (rdr.getCurrentLineNum() < 1)
            continue;
        /// FIXME: since vep 111, sv's insertion may become a sequences, e.g. tagttga...
        _varclass = rdr.getToken(columns["VARIANT_CLASS"]);
        this->lowercase(_varclass);
        if (isValidDNASequence(_varclass)) _varclass =  "insertion"; // seq -> category
        // if varType is "sv", "all", write all
        if ((_varType == "snv" || _varType == "indel") && (_varclass.find(_varType) == std::string::npos)) //
            continue;
        /// if the record is no in the queried input strains, skip.
        /// this ensure that we could identify impactful variants only from the queried strains.
        if (hasIND && (strainAbbrevs.size() > 0)) // if no IND column, then we need to parse the record
        {
            if (strainAbbrevs.hasIndex(rdr.getToken(columns["IND"])) < 0) continue;
        }
        /// if a heterozyote is annotated, treat them as the same to the reference genome, skip here
        /// FIXME: if don't exist ?
        if (rdr.getToken(columns["ZYG"]) == "HET")
            continue;
        
        /// aggregate results groupby location and transcript
        /// FIXME: what if transcript_id is '-'
        transcript_id = rdr.getToken(columns["Feature"]); 
        location = rdr.getToken(columns["Location"]);
        size_t tok0 = location.find_first_not_of("chr");
        location = location.substr(tok0, location.size());
        key = this->set_key(location, _varclass, _varType);
        // std::cout<<"Current keys: "<<key<<" <-->  "<<transcript_id<<std::endl;
        // std::cout<<rdr.getCurrentLineString()<<"\n\n"<<std::endl;

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
                std::string _csq = csq.substr(tokstart, dpos - tokstart);
                if (dpos > csz) dpos = csz;
                pVEP->consequence.addElementIfNew(_csq);
                tokstart = dpos + 1; // don't care if it overflows
            }
        }
        else
        {
            VEPSummary *pRecord = new VEPSummary(rdr.getToken(columns["Uploaded_variation"]), // alread trim # when read header
                                                 location, // rdr.getToken(columns["Location"]),
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
        if (codon == "-") 
        {
            // _expr = "";
            continue;
        }
        // std::string _codon(codon);
        upcase(codon);
        size_t pos = codon.find_first_of("/", 0);
        if (pos == std::string::npos) continue;
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
        Dynum<std::string> csq; // aggregate annotation to snp level, easy to remove duplicates
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
                            csqos <<"!";
                    }
                    csq.addElementIfNew(s);
                }
            } 
        }
        csqos << std::endl;
    }
    csqos.close();
}

void VarirantEeffectPredictor::writeVEPCsq(char* outFileName, bool prioritize) 
{
    std::ofstream csqos(outFileName);
    for (auto & k: this->keys)
    { // iter variant
        csqos << k;
        std::unordered_map<std::string, VEPSummary *>::iterator transxit = geneCodingMap[k].begin();
        /// aggreate snp annotation to gene level 
        std::unordered_map<std::string, CsqCoding> gene_csq; // used to store priority value
        for (; transxit != geneCodingMap[k].end(); transxit++)
        { // iter transcripts
            VEPSummary *pRecord = transxit->second;
            /// write records even symbol is "-" ? it's useless
            if (pRecord->symbol == "-") continue;  
            // aggregate annotation to gene level
            for (int i = 0; i < pRecord->consequence.size(); i++)
            {
                // std::string s = pRecord->symbol + "\t";
                std::string csq_str = pRecord->consequence.eltOf(i);
                // std::string csq_str2 = csq_str; // copy constructor
                // check
                if (CSQs.find(csq_str) == CSQs.end())
                {
                    std::cerr<< k << " " << pRecord->symbol <<"  ==> ";
                    std::cerr<< csq_str <<" <== is not recongnized. Skip."<<std::endl;
                    continue;
                }
                // current 
                CsqCoding csqcoding_cur = CsqCoding(pRecord->symbol, csq_str, CSQs[csq_str]);
                // update condon score
                if (csq_str.find("stop") != std::string::npos || 
                    csq_str == "missense_variant" || 
                    csq_str == "synonymous_variant" )
                {
                    csqcoding_cur.repr = this->codonChange(pRecord);
                } 
                // update to init
                csqcoding_cur.csqs.addElementIfNew(csqcoding_cur.repr);
                // aggregate results
                if ((gene_csq.find(pRecord->symbol) != gene_csq.end()))
                {
                    if (prioritize)
                    {  // max aggregate by impact score
                       if (csqcoding_cur.code > gene_csq[pRecord->symbol].code)
                            gene_csq[pRecord->symbol] = csqcoding_cur;
                    } 
                    else 
                    {  // do not prioritize
                        gene_csq[pRecord->symbol].raw.append("!"+csqcoding_cur.raw); // just concat strings
                        // gene_csq[pRecord->symbol].repr.append("!"+csqcoding_cur.repr);
                        gene_csq[pRecord->symbol].csqs.addElementIfNew(csqcoding_cur.repr); // remove duplicates
                        gene_csq[pRecord->symbol].code = std::max(gene_csq[pRecord->symbol].code, csqcoding_cur.code);
                    }
                }
                else 
                {  // init
                    gene_csq[pRecord->symbol] = csqcoding_cur;
                }
            } 
                   
        }
        // write gene level annotation 
        for (auto it = gene_csq.begin(); it != gene_csq.end(); ++it) 
        {
            // // csqos << "\t" << it->first << "\t" << it->second.raw;
            // write annotation
            csqos<< "\t" << it->first << "\t" << it->second.csqs.eltOf(0);
            for (int j = 1; j < it->second.csqs.size(); j++)
                csqos <<"!"<<it->second.csqs.eltOf(j); 
        }
        csqos << std::endl;
    }
    csqos.close();

}