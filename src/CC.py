
import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib
from pandas import Series

#Class component:

def class_component(quantitative_data_filename, data_process_origin, canopus_npc_df, filename_header, species_column, genus_column,family_column, metadata_df, reduced_df, min_class_confidence, min_recurrence, CC_component):
        
    if CC_component == True:
        row_ID_header = 'row ID'

        #read the file 

        quant_df_norm = pd.read_csv(quantitative_data_filename, sep=',')#,  index_col='row ID')
        quant_df_norm.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
        quant_df_norm.drop(list(quant_df_norm.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
        quant_df_norm.sort_index(axis=1, inplace=True)

        if data_process_origin == 'MZMine3':

            quant_df_norm.drop('row m/z', axis=1, inplace=True)
            quant_df_norm.drop('row retention time', axis=1, inplace=True)

            quant_df_norm.drop(['row ion mobility', 'correlation group ID', 'best ion', 'row ion mobility unit', 'row CCS', 
            'annotation network number', 'auto MS2 verify', 'identified by n=', 'partners', 'neutral M mass'], axis=1, inplace=True)
            quant_df_norm.set_index('row ID', inplace=True)

        else:

            quant_df_norm.drop('row m/z', axis=1, inplace=True)
            quant_df_norm.drop('row retention time', axis=1, inplace=True)
            quant_df_norm.set_index('row ID', inplace=True)

        quant_df_norm = quant_df_norm.apply(lambda x: x/x.max(), axis=0)
        
        quant_df_norm = quant_df_norm.transpose()
        quant_df_norm.reset_index(inplace=True)
        quant_df_norm.rename(columns={'index': filename_header}, inplace=True)
        quant_df_norm.set_index(filename_header, inplace=True)

        #retrieve the top 1 filename for each feature based on the quantitative table:
        df = quant_df_norm
        df = df.astype(float)
        df = df.apply(lambda s: s.abs().nlargest(1).index.tolist(), axis=1)
        df = df.to_frame()
        df[['row ID']] = pd.DataFrame(df[0].values.tolist(),index= df.index)
        df = df.drop([0], axis=1)
        #df.reset_index(inplace=True)
        df['row ID'] = pd.to_numeric(df['row ID'],errors = 'coerce')
        df['row ID'] = df['row ID'].fillna(0).astype(int)
        df.reset_index(inplace=True)

        #merged with the information from Canopus
        canopus_npc_df['shared name'] = canopus_npc_df['shared name'].astype(int) 
        df = pd.merge(df, canopus_npc_df[['shared name', 'NPC#class', 'NPC#class Probability']], how='left', left_on='row ID', right_on='shared name')#.dropna()
        df.rename(columns={'NPC#class Probability': 'NPC_class_Probability'}, inplace=True)
        df.drop('shared name', axis=1, inplace=True)

        #filter based on min_class_probability
        df = df[df.NPC_class_Probability > min_class_confidence]

        #calculare recurrence of each class 
        df= df[[filename_header, 'NPC#class']].groupby([filename_header,'NPC#class']).size().reset_index()
        df.rename(columns={0: 'recurrence'}, inplace=True)

        #filter based on minimun recurrence 
        df = df[df['recurrence'] >= min_recurrence].groupby(filename_header).agg(set)
        df.drop ('recurrence', axis=1, inplace=True)

        #add the species for each filename    
        df = pd.merge(metadata_df[[filename_header, species_column]], df, how='left', left_on=filename_header, right_on=filename_header).dropna()
        
        ## Retrieve classes reported in Lotus
        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv',sep=',').dropna()
        

        #create a set of species present in the metatada and reduce the lotus DB to it
        set_sp = set(metadata_df[species_column].dropna())
        LotusDB= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]

        #retrieve the chemical classes associated to the species and genus
        ccs = LotusDB[['organism_taxonomy_09species', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_09species').agg(set)
        ccs.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_species'}, inplace=True)
        ccg = LotusDB[['organism_taxonomy_08genus', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_08genus').agg(set)
        ccg.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_genus'}, inplace=True)

        #merge into a single dataframe
        cc = pd.merge(metadata_df[[filename_header, species_column, genus_column, family_column]], ccs, left_on= species_column, right_on='organism_taxonomy_09species', how='left')
        cc = pd.merge(cc,ccg, left_on= genus_column, right_on='organism_taxonomy_08genus', how='left')

        #obtain the difference between the predicted and reported compounds
        #first merge both dataframes
        df = pd.merge(df[[filename_header, 'NPC#class']],cc,on=filename_header, how='left')#.dropna()

        #check if the chemical classes from Sirius are reported in the species
        df['New_CC_in_sp'] = df["NPC#class"] - df["Chemical_class_reported_in_species"]  
        df['New_CC_in_genus'] = df["New_CC_in_sp"] - df["Chemical_class_reported_in_genus"]

        #function to check if there are new cc in ccs/ccg
        def is_empty(df):
                    """ function to check if the sets are empty or not 
                        Args:
                            df = Class component column CC 
                            Returns:
                            None
                    """
                    if df:
                        return 0.5 # if the column is not empty then 1, something is new in the sp &/ genus
                    else:
                        return 0

        df['CCs'] = df['New_CC_in_sp'].notnull().apply(is_empty)
        df['CCg'] = df['New_CC_in_genus'].notnull().apply(is_empty)
        df['CC'] = df['CCs'] + df['CCg']
        df['CC'] = df['CC'].fillna(1)

        #change all the NaN with a string to indicate lack of reports in the literature 
        df['Chemical_class_reported_in_species'] = df['Chemical_class_reported_in_species'].fillna('nothing in DB')
        df['Chemical_class_reported_in_genus'] = df['Chemical_class_reported_in_genus'].fillna('nothing in DB')
        df['New_CC_in_sp'] = df['New_CC_in_sp'].fillna('nothing in DB')
        df['New_CC_in_genus'] = df['New_CC_in_genus'].fillna('nothing in DB')

        df.to_csv('../data_out/CC_results.tsv', sep='\t')
        return df
    else:
        print ('No search was done because the Class component is not going to be calculated')

def class_component_ind_files(CC_component, repository_path, canopus_sample_suffix, metadata_df, filename_header, species_column, genus_column, family_column):
    """
    Function to recover the chemical classes predicted and reported from individual files and LOTUS accordingly, used for calculation of inventa non aligned data
    """
    if CC_component == True:
        df1 = pd.DataFrame()
        for r, d, f in os.walk(repository_path):
            for file in (f for f in f if f.endswith(canopus_sample_suffix)):
                complete_file_path =r+'/'+file
                read_file = pd.read_csv(complete_file_path, sep = '\t').dropna()
                read_file['partial_filename'] = read_file['name'].apply(lambda r: '_'.join(r.split('_')[:-3]) if len(r.split('_')) > 1 else r)
                read_file['partial_filename'] = read_file['partial_filename'].apply(lambda r: '_'.join(r.split('_')[1:]) if len(r.split('_')) > 1 else r)
                read_file = read_file[['partial_filename', 'class']].groupby('partial_filename').agg(set)
                df1 = df1.append(read_file, ignore_index=False)
            else:
                continue
        df1.reset_index(inplace=True)
        df2 = metadata_df
        s = [df2[df2[filename_header].str.contains(x)].index[0] for x in df1['partial_filename']]
        df1 = df1.assign(subcode=Series(data=df2[filename_header], index=s)).merge(df2[[filename_header, species_column]], left_on='subcode', right_on=filename_header).drop('subcode', axis='columns')
        df1.drop('partial_filename', axis=1, inplace=True)

        #now, try to get the reported chemical classes for the set 
        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv',sep=',').dropna()
        #create a set of species present in the metatada and reduce the lotus DB to it
        set_sp = set(df2[species_column].dropna())
        LotusDB= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]

        #retrieve the chemical classes associated to the species and genus
        df3 = LotusDB[['organism_taxonomy_09species', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_09species').agg(set)
        df3.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_species'}, inplace=True)
        df4 = LotusDB[['organism_taxonomy_08genus', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_08genus').agg(set)
        df4.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_genus'}, inplace=True)

        #merge into a single dataframe
        df5 = pd.merge(df2[[filename_header, species_column, genus_column, family_column]], df3, left_on= species_column, right_on='organism_taxonomy_09species', how='left')
        df5 = pd.merge(df5,df4, left_on= genus_column, right_on='organism_taxonomy_08genus', how='left')

        #obtain the difference between the predicted and reported compounds
        df = pd.merge(df1[[filename_header, 'class']],df5,on=filename_header, how='left').dropna()
        
        df['New_CC_in_sp'] = df["class"] - df["Chemical_class_reported_in_species"]  #check if the chemical classes from Sirius are reported in the species 
        df['New_CC_in_genus'] = df["New_CC_in_sp"] - df["Chemical_class_reported_in_genus"]
            
        #change all the NaN with a string to indicate lack of reports in the literature 
        df['Chemical_class_reported_in_species'] = df['Chemical_class_reported_in_species'].fillna('nothing in DB')
        df['Chemical_class_reported_in_genus'] = df['Chemical_class_reported_in_genus'].fillna('nothing in DB')

        df['New_CC_in_sp'] = df['New_CC_in_sp'].fillna('nothing in DB')
        df['New_CC_in_genus'] = df['New_CC_in_genus'].fillna('nothing in DB')

            
        def is_empty(df):
            """ function to check if the sets are empty or not 
                Args:
                    df = Class component column CC 
                    Returns:
                    None
            """
            if df:
                return 1 # if the column is not empty then 1, something is new in the sp &/ genus
            else:
                return 0
        df['CC'] = df['New_CC_in_sp'].apply(is_empty)
        df.to_csv('../data_out/LC_results.tsv', sep='\t')
        return df
    else:
        print ('No search was done because the Class component is not going to be calculated')