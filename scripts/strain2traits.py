#!/usr/bin/env python
# coding: utf-8
import sys, os, json
import requests
import numpy as np
import pandas as pd


class MPD:
    def __init__(self):
        """
        Vist here to get full API docs: https://phenome.jax.org/about/api#pheno_strainmeans
        
        See the data fields docs: https://phenome.jax.org/about/datafields#projects
        
        
        Example::
        ```
        mpd = MPD()
        projs = mpd.get_projects()
        # select interested data set    
        target_projs = projs.query('paneldesc == "inbred" & intervention == "nicotine")
        
        # data set example: from projsym field:
        projsym = "Vinyard1"
        # return a df (default) or json
        res = mpd.get_strain_means(projsym) 
        ```
        
        """
        self.host = "https://phenome.jax.org"
        
    def _get(self, url):
        """
        retrieve results
        """
        response = requests.get(url)
        if not response.ok:
            raise Exception("Retreve Error! Is input correct ?")          
        return response.json()
        
   
    def get_projects(self, custom_field=None, return_dataframe=True):
        """
        Fetch a list of MPD projects / data sets.
        """
        api = "/api/projects"
        if custom_field is not None:
            api += custom_field
        
        res = self._get(self.host + api)
        if return_dataframe: 
            return pd.DataFrame.from_dict(res['projects'])
        return res
    
    def get_projects_filter_by(self, filtername):
        """
        filtername: mpdsector, largecollab, panelsym, intervention
        
        """
        fn = str(filtername)
        assert fn.lower() in ['mpdsector', 'largecollab', 'panelsym', 'intervention']
        api = f"/api/project_filters/{filtername}"
        
        res = self._get(self.host + api)
        return res      
        
        
    def get_investigators(self, name=None):
        api = "/api/investigators"
        if name is not None:
            api += f"?name={name}"
            
        res = self._get(self.host + api)
        return res      
    
    def get_projsym_strains(self, projsym):
        """
        projsym:  project / data set identifier; usually the last name of first author followed by a digit
        
        """
        api = f"/api/projects/{projsym}/strains"
        res = self._get(self.host + api)
        return res  
    
    def get_projsym_publications(self, projsym):
        """
        projsym:  project / data set identifier; usually the last name of first author followed by a digit
        """
        api = f"/api/projects/{projsym}/strains"
    
        res = self._get(self.host + api)
        return res     
    
    def get_animaldata(self, measnum, covariate=None):
        """
        Fetch individual animal data for MPD strain survey phenotype measure(s) 
        identified by measnum (an integer or comma-separated list of multiple integers).
        
        covariate: the ID of another numeric measure involving the same animals
        """
        
        api = f"/api/pheno/animalvals/{measnum}"
        if covariate is not None:
            api = api + f"?covariate={covariate}"
        
        res = self._get(self.host + api)
        return res           
        
    def get_animaldata_series(self, measnum):
        """
        Fetch individual animal data for a series of measures identified by measnum.
        """
        api = f"/api/pheno/animalvals/series/{measnum}"
        res = self._get(self.host + api)
        return res    
    
    def get_strain_means(self, selector, return_dataframe=True):
        """
        Get unadjusted strain means for one or more MPD phenotype measures.
        
        selector:  
           - measnum: seperate each id by comma
           - projsym: a valid project symbol such as Vinyard1. 
        
        """
        api = f"/api/pheno/strainmeans/{selector}"
        res = self._get(self.host + api)
        if return_dataframe: 
            return pd.DataFrame.from_dict(res['strainmeans'])
        
        return res     
    
    def get_strain_lsmeans(self, selector, return_dataframe=True):
        """Get model-adjusted least-square strain means for one or more 
           MPD phenotype measures in the "lsmeans" result element. 
           
          selector:  
             - measnum: seperate each id by comma
             - projsym: a valid project symbol such as Vinyard1. 
           
           Least squares means are arithmetic means adjusted by model term(s). 
           They represent a better estimate of the true population mean than the unadjusted group means,
           and are less sensitive to missing data. 
           
           It may be noticed that LSM SEMs tend to be uniform when the model is balanced.
        
        """
        api = f"/api/pheno/lsmeans/{selector}"
        res = self._get(self.host + api)
        if return_dataframe: 
            return pd.DataFrame.from_dict(res['lsmeans'])
        
        return res      
    
    
    
    def get_measureinfo(self, selector):
        """
        Get descriptions, units, and other metadata for one or more MPD measures.
        
        selector:  
           - measnum: seperate each id by comma
           - projsym: a valid project symbol such as Vinyard1. 
        """
        api = f"/api/pheno/measureinfo/{selector}"
        res = self._get(self.host + api)
        return res   
    
    def get_seriesmembers(self, keymeasnum):
        """
        Get measure IDs, descriptions, and other metadata for all members of a measure series identified by keymeasnum.
        """
        api = f"/api/pheno/seriesmembers/{keymeasnum}"
        res = self._get(self.host + api)
        return res  
    
    def get_measures_by_ontology(self, ont_term):
        """
        The default use of this endpoint is to get measure IDs, descriptions,
        and other metadata for all measures that have been annotated to MP, VT, or MA ontology
        
        """
        api = f"/api/pheno/measures_by_ontology/{ont_term}"
        
        res = self._get(self.host + api)
        return res    


def strain2trait(trait, strain, measnum, outdir, use_rawdata=False):
    mpd = MPD()
    ## parse ids already    
    IDS_ORIG = measnum
    ## ###
    # required file: SNP database metadata 
    # should have these columns:
    # strain_fullname, strain_abbr, mpd_strainid
   
    if strain.endswith(".txt"):
        strains = pd.read_table(strain, dtype=str)
    else:
        strains = pd.read_csv(strain, dtype=str)
    
    # contain columns:
    # measnum, strain, strainid, mean

    # FIXME: makesure all measnums existed in MPD
    if (trait is None) or trait == "" or use_rawdata:
        if measnum.find("-") != -1: 
            measnum = measnum.split("-")[0]
        
        if use_rawdata:
            traits = mpd.get_animaldata(measnum)
        else:
            traits = mpd.get_strain_means(measnum)  
        # FIXME: Na could not convert to str
        traits = traits.astype(str)
    elif trait.endswith(".txt"):
        traits  = pd.read_table(trait, dtype=str)
    elif trait.endswith(".csv"):
        traits  = pd.read_csv(trait, dtype=str)
    else:
        raise Exception("Trait file Error")

    if 'measnumSex' in traits.columns:
        measumid = traits.columns.get_loc('measnumSex')
        IDS = IDS_ORIG
    elif 'measnum' in traits.columns:
        measumid = traits.columns.get_loc('measnum')
        IDS = measnum
    else:
        sys.exit("Missing column measnum or measnumSex!") 
        
    traits = traits[traits.iloc[:, measumid] == IDS]
    if traits.shape[0] == 0:
        sys.exit("ERROR! NOT Matched trait ids found !!!")
    
    mpd2strain = {v.loc['MPD'] : v.loc['Abbr'] for k,v in strains.iterrows()}
    # remove extra " ", it's no easy to detect this bug!
    traits['strainid'] = traits.strainid.str.strip()
    traits['strain_abbr'] = traits.strainid.map(mpd2strain) 
    traits = traits[traits.strain_abbr.isin(strains['Abbr'])]
    # drop strains could not map to SNP database
    traits.dropna(subset=['strain_abbr'], inplace=True)
    trait_ids = traits.iloc[:, measumid].unique()

    if len(trait_ids) == 0:
        sys.exit("ERROR! NOT Matched trait ids found !!!")
    if len(trait_ids) != 1:
        print(f"Warning! ID Not Matched! Check file: {IDS}", file=sys.stderr)

    ## split male and female
    def write(temp, meta, tid, field):
        os.makedirs(os.path.join(outdir,f"MPD_{tid}"), exist_ok=True)
        temp.loc[:,['strain_abbr','strain']].to_csv(os.path.join(outdir,f"MPD_{tid}/strain.{tid}.txt"),
                                                    index=False, header=None, sep="\t")
        if meta['measures_info'][0]['datatype'] == 'categorical':
            print("Cateogorical measure found!")
            catout = os.path.join(outdir,f"MPD_{tid}/trait.{tid}.categorical")
            os.system(f"touch {catout}")
            temp[field] = temp[field].astype('category')
            temp.loc[:,['strain_abbr', 'strain', field]].to_csv(os.path.join(outdir,f"MPD_{tid}/trait.{tid}.txt"),
                                                      index=False, header=None, sep="\t")
        else:
            temp.loc[:,['strain_abbr', 'strain', field]].to_csv(os.path.join(outdir,f"MPD_{tid}/trait.{tid}.txt"), 
                                                      index=False, header=None, sep="\t", float_format='%.7f')       
    data = 'value' if use_rawdata else 'mean'
    tid = IDS_ORIG
    temp = traits.copy()
    # get metadata from MPD
    if tid.find("-") != -1:
        m, s = tid.split("-") # measum and sex
        meta = mpd.get_measureinfo(selector=m)
        sex = temp.sex.unique()
        if s in sex:
            write(temp.query(f'sex == "{s}"'), meta, f"{m}-{s}", data)
        else:
            print(f"the suffix of {tid} represent sex? Should be: f or m ")
        return
    else:
        meta = mpd.get_measureinfo(selector=tid)
    #breakpoint()
    # split male and female
    if (tid.find("-") == -1) and (temp.sex.nunique() == 2):
        sex = temp.sex.unique()
        for s in sex: 
            os.makedirs(os.path.join(outdir,f"MPD_{tid}-{s}"), exist_ok=True)
            write(temp.query(f'sex == "{s}"'), meta, f"{tid}-{s}", data)

        if len(sex) == 1: return # don't write duplicate outputs

    ## if only one sex present, ignore splitted id.
    # temp[data] = temp[data].astype(float)
    # temp = temp.groupby(['strain_abbr', 'strain'], as_index=False)[data].mean()

    # both sex exits, write an aggregated output
    os.makedirs(os.path.join(outdir,f"MPD_{tid}"), exist_ok=True)
    write(temp, meta, tid, data)



strain2trait(snakemake.params['trait'], snakemake.input['strain'], 
             snakemake.params['traitid'], snakemake.params['outdir'], 
             snakemake.params['rawdata'])


# if __name__ == "__main__":
#     strain2trait("", strain="/data/bases/shared/haplomap/PELTZ_20200429/strains.metadata.csv",
#                  measnum= "31870",
#                  outdir="/data/bases/shared/haplomap/PELTZ_20200429/test")