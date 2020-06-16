#!/usr/bin/env python
# coding: utf-8
import sys, os
import numpy as np
import pandas as pd

def strain2trait(trait, strain, ids, outdir, categorical=True):

    # required file
    # should have these columns:
    # strain_fullname, strain_abbr, MPDid
    if strain.endswith(".txt"):
        strains = pd.read_table(strain, dtype=str)
    else:
        strains = pd.read_csv(strain, dtype=str)
    
    # contain columns:
    # measnum, strain, mpdid(strainid), mean
    if trait.endswith(".txt"):
        traits  = pd.read_table(trait, dtype=str)
    else:
        traits  = pd.read_csv(trait, dtype=str)
    ## parse ids already    
    IDS = ids
    ## 
    if 'measnumSex' in traits.columns:
        measumid = traits.columns.get_loc('measnumSex')
    elif 'measnum' in traits.columns:
        measumid = traits.columns.get_loc('measnum')
    else:
        sys.exit("Missing column measnum or measnumSex!") 
        
    traits = traits[traits.iloc[:, measumid] == IDS]
    if traits.shape[0] == 0:
        sys.exit("ERROR! NOT Matched trait ids found !!!")


    # FIXME: handle `strain` B10 , it does not have a MPD strainid!
    # use JAX id instead
    if pd.isna(strains.loc[strains.Abbr == 'B10','MPD']).any().sum() or strains.loc[strains.Abbr == 'B10','MPD'] == "":
        strains.loc[strains.Abbr == 'B10','MPD'] = '462'  
        strains['MPD'] = strains['MPD'].astype(float).astype(pd.Int16Dtype()).astype(str)
        strains.loc[strains.Abbr == 'B10','MPD'] = '000462' 
    # handle b10 strain id!
    traits.loc[ traits['strain'].isin(['B10.D2-H2<d>/n2SnJ','B10']), 'strainid'] = '000462'
    #
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
    if len(trait_ids) != len(np.unique(IDS)):
        print(f"Warning! Some IDs Not Matched! Skipped! Check file: {ids}", file=sys.stderr)

    for tid in trait_ids:
        temp = traits[traits.iloc[:,measumid] == tid]
        os.makedirs(os.path.join(outdir,f"MPD_{tid}"), exist_ok=True)

        temp.loc[:,['strain_abbr','strain']].to_csv(os.path.join(outdir,f"MPD_{tid}/strain.{tid}.txt"),
                                                    index=False, header=None, sep="\t")
        
        # FIXME: if categorical
        if categorical and temp['mean'].nunique() < temp['strain_abbr'].nunique():
            print("Cateogorical measure found!")
            catout = os.path.join(outdir,f"MPD_{tid}/trait.{tid}.categorical")
            os.system(f"touch {catout}")
            temp['mean'] = temp['mean'].astype('category')
            temp.loc[:,['strain_abbr','mean']].to_csv(os.path.join(outdir,f"MPD_{tid}/trait.{tid}.txt"),
                                                      index=False, header=None, sep="\t")
        else:
            temp.loc[:,['strain_abbr','mean']].to_csv(os.path.join(outdir,f"MPD_{tid}/trait.{tid}.txt"), 
                                                      index=False, header=None, sep="\t", float_format="%.7f")

strain2trait(snakemake.input['trait'], snakemake.input['strain'], 
             snakemake.params['traitid'], snakemake.params['outdir'], 
             snakemake.params['has_categorical'])