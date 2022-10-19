
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
from tqdm import tqdm

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
        df = quant_df_norm.transpose()
        df = df.astype(float)
        df = df.apply(lambda s: s.abs().nlargest(1).index.tolist(), axis=1)
        df = df.to_frame()
        df[['filename']] = pd.DataFrame(df[0].values.tolist(),index= df.index)
        df = df.drop([0], axis=1)
        #df.reset_index(inplace=True)
        #df['filename'] = pd.to_numeric(df['filen'],errors = 'coerce')
        #df['row ID'] = df['row ID'].fillna(0).astype(int)
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

#class component general for NPClassifier

def class_component_ind_files(CC_component, repository_path, ionization_mode, min_class_confidence, metadata_df, filename_header, species_column, genus_column, family_column):
    """
    Function to recover the chemical classes predicted and reported from individual files and LOTUS accordingly, used for calculation of inventa non aligned data
    """
    if CC_component == True:
        
        path = os.path.normpath(repository_path)
        samples_dir = [directory for directory in os.listdir(path)]

        df =pd.DataFrame()
        files = []
        for directory in tqdm(samples_dir):
            canopus_path = os.path.join(path, directory, ionization_mode, directory + '_WORKSPACE_SIRIUS', 'canopus_formula_summary_adducts.tsv')
            try:
                canopus_df =pd.read_csv(canopus_path, sep='\t')
            except FileNotFoundError:
                continue
            except NotADirectoryError:
                continue

            #recover filenames 
            files.append(directory)
            
            #read canpus filename
            canopus_df =pd.read_csv(canopus_path, sep='\t')
            canopus_df['partial_filename'] = canopus_df['id'].apply(lambda r: '/'.join(r.split('/')[:1]) if len(r.split('_')) > 1 else r)
            canopus_df['partial_filename'] = canopus_df['partial_filename'].apply(lambda r: '_'.join(r.split('_')[:-3]) if len(r.split('_')) > 1 else r)
            canopus_df['partial_filename'] = canopus_df['partial_filename'].apply(lambda r: '_'.join(r.split('_')[1:]) if len(r.split('_')) > 1 else r)

            canopus_df = canopus_df[['partial_filename','id', 'molecularFormula', 'adduct', 'NPC#pathway',
                'NPC#pathway Probability', 'NPC#superclass',
                'NPC#superclass Probability', 'NPC#class', 'NPC#class Probability']]
            canopus_df.rename(columns={'NPC#class Probability': 'NPC_class_Probability'}, inplace=True) 
            canopus_df['shared name'] = canopus_df['id'].str.split('_').str[-1].astype(int)
            canopus_df.drop('id', axis=1, inplace=True)
            canopus_df.rename(columns={'shared name': 'row ID', 'adduct': 'adduct (sirius)', 'molecularFormula': 'MF (sirius)', 'name': 'Compound name (sirius)'}, inplace=True) 
            canopus_df.drop(canopus_df[canopus_df.NPC_class_Probability > min_class_confidence].index, inplace=True)
            canopus_df.drop(['NPC_class_Probability', 'NPC#superclass Probability', 'NPC#pathway Probability'], axis=1, inplace=True)
            canopus_df = canopus_df[['partial_filename', 'NPC#class']].groupby('partial_filename').agg(set)
            df = df.append(canopus_df, ignore_index=False)

        df.reset_index(inplace=True)
        df = pd.merge(df, metadata_df[[filename_header, species_column]], how= 'left', left_on='partial_filename', right_on=filename_header)
        df.drop('partial_filename', axis=1, inplace=True)

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
        
        #if we considered lack of reports in DB... all the proposed cc will be new to that particular species, hence.. the real value of CC is 1, not zero
        string = 'nothing in DB'
        df.loc[df['Chemical_class_reported_in_species'] ==string, 'CC'] = 1
       
        pathout = os.path.join(path, 'results/')
        os.makedirs(pathout, exist_ok=True)
        pathout = os.path.join(pathout, 'Class_component_results' +'_' + ionization_mode + '.tsv')
        df.to_csv(pathout, sep ='\t')
        return df
    else:
        print ('No search was done because the Class component is not going to be calculated')


#class component general for NPClassifier PF1600

def class_component_ind_files_PF1600(CC_component, repository_path, ionization_mode, min_class_confidence, metadata_df, filename_header, species_column, genus_column, family_column):
    """
    Function to recover the chemical classes predicted and reported from individual files and LOTUS accordingly, used for calculation of inventa non aligned data
    """
    if CC_component == True:
        
        path = os.path.normpath(repository_path)
        samples_dir = [directory for directory in os.listdir(path)]

        df =pd.DataFrame()
        files = []
        for directory in tqdm(samples_dir):
            canopus_path = os.path.join(path, directory, ionization_mode, directory + '_WORKSPACE_SIRIUS', 'npc_summary.csv')
            #canopus_path = os.path.join(path, directory, ionization_mode, directory + '_WORKSPACE_SIRIUS', 'canopus_summary_adducts.tsv')
            try:
                canopus_df =pd.read_csv(canopus_path, sep=',')
            except FileNotFoundError:
                continue
            except NotADirectoryError:
                continue

            #recover filenames 
            files.append(directory)
            
            #read canpus filename
            canopus_df =pd.read_csv(canopus_path, sep=',')
            canopus_df['partial_filename'] = canopus_df['directoryName'].apply(lambda r: '/'.join(r.split('/')[8:]) if len(r.split('_')) > 1 else r)
            canopus_df['partial_filename'] = canopus_df['partial_filename'].apply(lambda r: '_'.join(r.split('_')[:-3]) if len(r.split('_')) > 1 else r)
            canopus_df['partial_filename'] = canopus_df['partial_filename'].apply(lambda r: '_'.join(r.split('_')[1:]) if len(r.split('_')) > 1 else r)
            canopus_df.rename(columns={'name': 'shared name','classProbability': 'NPC_class_Probability', 'class':'NPC#class'}, inplace=True)

            canopus_df = canopus_df[['partial_filename', 'shared name', 'NPC#class', 'NPC_class_Probability']]
            canopus_df = canopus_df[canopus_df.NPC_class_Probability > min_class_confidence]
            canopus_df = canopus_df[['partial_filename', 'NPC#class']].groupby('partial_filename').agg(set)
            df = df.append(canopus_df, ignore_index=False)

        df.reset_index(inplace=True)
        df = pd.merge(df, metadata_df[[filename_header, species_column]], how= 'left', left_on='partial_filename', right_on=filename_header)
        df.drop('partial_filename', axis=1, inplace=True)

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
        
        #if we considered lack of reports in DB... all the proposed cc will be new to that particular species, hence.. the real value of CC is 1, not zero
        string = 'nothing in DB'
        df.loc[df['Chemical_class_reported_in_species'] ==string, 'CC'] = 1
       
        pathout = os.path.join(path, 'results/')
        os.makedirs(pathout, exist_ok=True)
        pathout = os.path.join(pathout, 'Class_component_results' +'_' + ionization_mode + '.tsv')
        df.to_csv(pathout, sep ='\t')
        return df
    else:
        print ('No search was done because the Class component is not going to be calculated')