
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

def class_component(df3, filename_header, species_column, genus_column,family_column, metadata_df, reduced_df, col_id_unique, min_class_confidence, min_recurrence, CC_component):
    """ function to compute the class component based on the possible presence of new chemical classes 
    Args:
        df3: canopus_sunmmary
        Returns:
        None
    """

    if CC_component == True:

        def top_ions(col_id_unique, reduced_df):

                """ function to compute the top species, top filename and top species/plant part for each ion 
            Args:
                df1 = reduced_df, table of with index on sp/part column and features only.
                df2 = quantitative.csv file, output from MZmine
                Returns:
                None
            """
            #computes the % for each feature
                dfA = reduced_df.copy()
                #dfA = dfA.transpose()
                dfA = dfA.div(dfA.sum(axis=1), axis=0)
                dfA.reset_index(inplace=True)
                dfA.rename(columns={'index': 'row ID'}, inplace=True)
                dfA.set_index('row ID', inplace=True)
                dfA = dfA.astype(float)
                dfA['Feature_specificity'] = dfA.apply(lambda s: s.abs().nlargest(1).sum(), axis=1)
                dfA.reset_index(inplace=True)
                #df1 = df1.drop([0], axis=1)
                dfA = dfA[['row ID', 'Feature_specificity']]
                dfA['row ID']=dfA['row ID'].astype(int)

                #computes the top filename for each ion 
                df2 = reduced_df.copy()
                df2 = df2.div(df2.sum(axis=1), axis=0)
                df2 = df2.copy()
                df2 = df2.astype(float)
                df2 = df2.apply(lambda s: s.abs().nlargest(1).index.tolist(), axis=1)
                df2 = df2.to_frame()
                df2['filename'] = pd.DataFrame(df2[0].values.tolist(), index= df2.index)
                df2 = df2.drop([0], axis=1)
                df2.reset_index(inplace=True)
                df2.rename(columns={'index': 'row ID'}, inplace=True)

                df = pd.merge(left=dfA,right=df2, how='left',on='row ID')

                if col_id_unique != 'filename':
                    #computes the top species/part for each feature 
                    df3 = reduced_df.copy()
                    df3 = df3.transpose()
                    df3 = df3.astype(float)
                    df3 = df3.apply(lambda s: s.abs().nlargest(1).index.tolist(), axis=1)
                    df3 = df3.to_frame()
                    df3[[col_id_unique]] = pd.DataFrame(df3[0].values.tolist(),index= df3.index)
                    df3 = df3.drop([0], axis=1)
                    df3.reset_index(inplace=True)
                    df3.rename(columns={'index': 'row ID'}, inplace=True)
                    df3['row ID'] = df3['row ID'].astype(int)
                    
                    #merge all the data 
                    df = pd.merge(left=df3, right=df, how='left', on='row ID')
                else: 
                    df
                return df
                
        df1 = top_ions(col_id_unique, reduced_df)
        df2 = metadata_df.copy()

        df3 = get_canopus_pred_classes(canopus_npc_summary_filename, CC_component)

        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv',sep=',').dropna()

        # merge with top filename with iones 
        #df3['shared name'] = df3['name'].str.split('_').str[-1].astype(int) #use this line if you don't have a column with the name/shared name

        ## Retrieve classes predicted by Sirius         
        #the specificity_df is used to assign the main biological source to each feature. 
        df3.rename(columns={'name': 'shared name'}, inplace=True)
        df3 = pd.merge(left=df1[['row ID', filename_header]], right=df3[['shared name', 'class', 'classProbability']], how='left', left_on='row ID', right_on='shared name').dropna()
        df3.drop('shared name', axis=1, inplace=True)

        #filter based on min_class_probability
        df3.drop(df3[df3.classProbability > min_class_confidence].index, inplace=True)

        #calculare recurrence of each class 
        df4= df3[[filename_header, 'class']].groupby([filename_header,'class']).size().reset_index()
        df4.rename(columns={0: 'recurrence'}, inplace=True)

        df4 = df4[df4['recurrence'] >= min_recurrence].groupby(filename_header).agg(set)
        df4.drop ('recurrence', axis=1, inplace=True)
            
        df5 = pd.merge(left=df2[[filename_header, species_column]], right=df4, how='left', left_on=filename_header, right_on=filename_header).dropna()

        ## Retrieve classes reported in Lotus

        #create a set of species present in the metatada and reduce the lotus DB to it
        set_sp = set(df5[species_column].dropna())
        LotusDB= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]

        #retrieve the chemical classes associated to the species and genus
        df7 = LotusDB[['organism_taxonomy_09species', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_09species').agg(set)
        df7.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_species'}, inplace=True)
        df8 = LotusDB[['organism_taxonomy_08genus', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_08genus').agg(set)
        df8.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_genus'}, inplace=True)
        #merge into a single dataframe
        df9 = pd.merge(df2[[filename_header, species_column, genus_column, family_column]], df7, left_on= species_column, right_on='organism_taxonomy_09species', how='left')
        df9 = pd.merge(df9,df8, left_on= genus_column, right_on='organism_taxonomy_08genus', how='left')
        #df9.head(5)
        #df9.to_csv('../data_out/reported_classes_df.tsv', sep='\t')

        #obtain the difference between the predicted and reported compounds
        df = pd.merge(df5[[filename_header, 'class']],df9,on=filename_header, how='left').dropna()

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
        df['CC'] = df['CC'].fillna(1)

        #df = pd.merge(df2[[filename_header]], df,how= 'left', on=filename_header)
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