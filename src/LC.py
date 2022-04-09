import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib

#literature component
 
def literature_component(LC_component, max_comp_reported_sp, max_comp_reported_g, max_comp_reported_f):
    """ function to compute the literature component based on the metadata and combinend information of the Dictionary of natural products and the Lotus DB, 
    Args:
        df2 = metadata_df

    Returns:
        None
    """
    if LC_component == True:
        df = pd.read_csv('../data_out/metadata_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv', 
                       sep=',').dropna()

        #create a set of species from the metadata table
        set_sp = set(df['ATTRIBUTE_Species'].dropna()) #dropna is used to erase all the QC, blanks, etc not having a species associated
        #create a set of genera from the metadata table 
        set_g = set(df['ATTRIBUTE_Genus'].dropna()) #dropna is used to erase all the QC, blanks, etc not having a species associated
        #create a set of families from the metadata table 
        set_f = set(df['ATTRIBUTE_Family'].dropna()) #dropna is used to erase all the QC, blanks, etc not having a species associated

        #reduce LotusDB to sp in set 
        LotusDBs= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]
        LotusDBs = LotusDBs[['organism_name','organism_taxonomy_09species','Reported_comp_Species']].drop_duplicates()
        LotusDBs['Reported_comp_Species'] =LotusDBs['Reported_comp_Species'].astype(int)

        #reduce LotusDB to genera in set 
        LotusDBg= LotusDB[LotusDB['organism_taxonomy_08genus'].isin(set_g)]
        LotusDBg = LotusDBg[['organism_taxonomy_08genus','Reported_comp_Genus']].drop_duplicates()
        LotusDBg['Reported_comp_Genus'] =LotusDBg['Reported_comp_Genus'].astype(int)

        #reduce LotusDB to families in set 
        LotusDBf= LotusDB[LotusDB['organism_taxonomy_06family'].isin(set_f)]
        LotusDBf = LotusDBf[['organism_taxonomy_06family','Reported_comp_Family']].drop_duplicates()
        LotusDBf['Reported_comp_Family'] =LotusDBf['Reported_comp_Family'].astype(int)
        #LotusDB.head()
        
        df = pd.merge(df[['filename', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species']],
                LotusDBs[['organism_taxonomy_09species', 'Reported_comp_Species']],
                how= 'left', left_on='ATTRIBUTE_Species', right_on='organism_taxonomy_09species').drop_duplicates(subset=['filename'])
        df.drop('organism_taxonomy_09species', axis=1, inplace=True)
        df = pd.merge(df,
                LotusDBg,
                how= 'left', left_on='ATTRIBUTE_Genus', right_on='organism_taxonomy_08genus').drop_duplicates(subset=['filename'])
        df.drop('organism_taxonomy_08genus', axis=1, inplace=True)
        df = pd.merge(df,
                LotusDBf,
                how= 'left', left_on='ATTRIBUTE_Family', right_on='organism_taxonomy_06family').drop_duplicates(subset=['filename'])
        df.drop('organism_taxonomy_06family', axis=1, inplace=True)


        df = df.fillna(0) #assumign species not present in LotusDB the number of reported compounds is set to 0
        df['LC'] = 1-df['Reported_comp_Species'].div(max_comp_reported_sp*100)-df['Reported_comp_Genus'].div(max_comp_reported_g*100)-df['Reported_comp_Family'].div(max_comp_reported_f*100)
        #df['LC'] = df.apply(literature_report, axis=1)
        df.to_csv('../data_out/LC_results.tsv', sep='\t')
        return df

    else:
        print('Literature component not calculated')
        