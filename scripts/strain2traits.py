#!/usr/bin/env python
# coding: utf-8
import sys, os
import numpy as np
import pandas as pd

def strain2trait(trait, strain, ids, outdir):
    strains = pd.read_csv(strain)
    traits  = pd.read_csv(trait)
    IDS = ids
    # measum ids colname
    measumid = traits.columns.get_loc('measnumSex')
    traits = traits[traits.iloc[:,measumid] == IDS]
    if traits.shape[0] == 0:
        sys.exit("ERROR! NOT Matched trait ids found !!!")

    mpd2strain = {v.loc['MPD'] : v.loc['Abbr'] for k,v in strains.iterrows()}
    traits['strain_abbr'] = traits.strainid.map(mpd2strain) 
    traits = traits[traits.strain_abbr.isin(strains['Abbr'])]
    # drop strains could not map to SNP database
    traits.dropna(subset=['strain_abbr'], inplace=True)
    trait_ids = traits.iloc[:,measumid].unique()
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
        if temp['mean'].unique().shape[0] < temp['strain_abbr'].unique().shape[0]:
            print("Cateogorical measure found!")
            catout = os.path.join(outdir,f"MPD_{tid}/strain.{tid}.categorical")
            os.system(f"touch {catout}")
            temp['mean'] = temp['mean'].astype(int)
            temp.loc[:,['strain_abbr','mean']].to_csv(os.path.join(outdir,f"MPD_{tid}/trait.{tid}.txt"),
                                                      index=False, header=None, sep="\t")
        else:
            temp.loc[:,['strain_abbr','mean']].to_csv(os.path.join(outdir,f"MPD_{tid}/trait.{tid}.txt"), 
                                                      index=False, header=None, sep="\t", float_format="%.7f")

strain2trait(snakemake.input['trait'], snakemake.input['strain'], 
             snakemake.params['traitid'],   snakemake.params['outdir'])