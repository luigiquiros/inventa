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
 
def literature_component(df, LC_component, min_comp_reported, max_comp_reported):
    """ function to compute the literature component based on the metadata and combinend information of the Dictionary of natural products and the Lotus DB, 
    Args:
        df2 = metadata_df

    Returns:
        None
    """
    if LC_component == True:
        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv', 
                       sep=',').dropna()

        #create a set of species from the metadata table
        set_sp = set(df['ATTRIBUTE_Species'].dropna()) #dropna is used to erase all the QC, blanks, etc not having a species associated

        #reduce LotusDB a sp in set 
        LotusDB= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]
        LotusDB = LotusDB[['organism_name',
            'organism_taxonomy_06family', 'organism_taxonomy_08genus',
            'organism_taxonomy_09species', 'Reported_comp_Family',
            'Reported_comp_Genus', 'Reported_comp_Species']].drop_duplicates()
        #LotusDB.head()
        
        df = pd.merge(df[['filename', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species']],
                LotusDB[['organism_taxonomy_09species', 'Reported_comp_Family','Reported_comp_Genus', 'Reported_comp_Species']],
                how= 'left', left_on='ATTRIBUTE_Species', right_on='organism_taxonomy_09species').drop_duplicates(subset=['filename'])
        df.drop('organism_taxonomy_09species', axis=1, inplace=True)
        df = df.fillna(0) #assumign species not present in LotusDB the number of reported compounds is set to 0
        df['Reported_comp_Species'] = df['Reported_comp_Species'].astype(int) 

        
        def literature_report(df):
            """ function to give a weigth to the counts of the reported compouds according to the used
            Args:
                df = Literature_component output
            Returns:
                None
            """
            if (df['Reported_comp_Species'] <= min_comp_reported):
                return 1
            elif (df['Reported_comp_Species'] <= max_comp_reported & df['Reported_comp_Species'] >= min_comp_reported): 
                return 0.5
            else:
                return 0

        df['LC'] = df.apply(literature_report, axis=1)
        df.to_csv('../data_out/LC_results.tsv', sep='\t')
        return df

    else:
        print('Literature component not calculated')
        