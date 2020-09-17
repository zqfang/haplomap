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


    def add_suffix(self, measnum):

        # contain columns:
        # measnum, strain, strainid, mean
        outids = []
        sex = self.get_measureinfo(measnum)['measures_info'][0]['sextested']
        if sex == "both":
            outids += [measnum+"-f", measnum+"-m"]
        elif sex == "m":
            outids.append(measnum+"-m")
        elif sex == "f":
            outids.append(measnum+"-f")
        else:
            print(f"MPD id: {measnum}, sex field is {sex}")
        return outids
        


if __name__ == "__main__":
        ## ###
   
    ## usage: generate_id.py measurm.txt neasum_suffixed.txt
    traits = sys.argv[1]
    outfile = sys.argv[2]
    with open(traits,'r') as t:
        trait_ids = t.read().strip().split()
    mpd = MPD()
    newids = []
    for idx in trait_ids:
        newids += mpd.add_suffix(idx)
    
    with open(outfile,'w') as out:
        for n in newids:
            out.write(n+"\n")

