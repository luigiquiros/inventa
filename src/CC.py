
import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib

#Class component:

def sirius_classes1(df1,df2,df3, min_recurrence, CC_component): 
    """ function to find the chemical classes proposed by sirius and assign them to a specific species based on the top specificity of the feature
    Args:
        df1 = specificity_df
        df2 = metadata_df
        df3 = output from SIRIUS + Canopus 

    Returns:
        None
    """
    if CC_component == True:
        # merge with top filename with iones 
        #df3['shared name'] = df3['name'].str.split('_').str[-1].astype(int) #use this line if you don't have a column with the name/shared name
        
        #the specificity_df is used to assign the main biological source to each feature. 
        df3.rename(columns={'name': 'shared name'}, inplace=True)
        df3 = pd.merge(left=df1[['row ID', 'filename']], right=df3[['shared name', 'class', 'classProbability']], how='left', left_on='row ID', right_on='shared name').dropna()
        df3.drop('shared name', axis=1, inplace=True)

        #filter based on min_class_probability
        df3.drop(df3[df3.classProbability < min_class_confidence].index, inplace=True)

        #calculare recurrence of each class 
        df4= df3[['filename', 'class']].groupby(['filename','class']).size().reset_index()
        df4.rename(columns={0: 'recurrence'}, inplace=True)

        df4 = df4[df4['recurrence'] >= min_recurrence].groupby('filename').agg(set)
        df4.drop ('recurrence', axis=1, inplace=True)
        
        df = pd.merge(left=df2[['filename', 'ATTRIBUTE_Species']], right=df4, how='left', left_on='filename', right_on='filename').dropna()
        df.drop('ATTRIBUTE_Species', axis=1, inplace=True)
        return df
    else:
        print ('No search was done because the Class component is not going to be calculated')

def search_reported_class(df, CC_component):
    """ function to search the reported chemical classes in each species of the set 
    Args:
        df = metadata_df
        Returns:
        None
    """
    if CC_component == True:
        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv',
                        sep=',').dropna()
        
        #create a set of species present in the metatada and reduce the lotus DB to it
        set_sp = set(df['ATTRIBUTE_Species'].dropna())
        LotusDB= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]

        #retrieve the chemical classes associated to the species and genus
        df4 = LotusDB[['organism_taxonomy_09species', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_09species').agg(set)
        df4.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_species'}, inplace=True)
        df5 = LotusDB[['organism_taxonomy_08genus', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_08genus').agg(set)
        df5.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_genus'}, inplace=True)

        #merge into a single dataframe

        df = pd.merge(df[['filename', 'ATTRIBUTE_Species', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Family', 'ATTRIBUTE_Family']],df4,left_on= 'ATTRIBUTE_Species', right_on='organism_taxonomy_09species', how='left')
        df = pd.merge(df,df5,left_on= 'ATTRIBUTE_Genus', right_on='organism_taxonomy_08genus', how='left') 
        return df
    else:
        print ('No search was done because the Class component is not going to be calculated')

def class_component(df1, df2, df3, CC_component):
    """ function to compute the class component based on the possible presence of new chemical classes 
    Args:
        df1 = reported_classes_df 
        df2 = sirius_classes_df
        df3 = metadata_df
        Returns:
        None
    """
    if CC_component == True:
    #merge the both tables
        df = pd.merge(df1,df2,on='filename', how='left').dropna()

        #get the difference between sets 

        df['New_in_species'] = df["class"] - df["Chemical_class_reported_in_species"]  #check if the chemical classes from Sirius are reported in the species
        df['New_in_genus'] = df["New_in_species"] - df["Chemical_class_reported_in_genus"]  #check if the NEW chemical classes in the species are reported in the genus

        #Add the weight accordingly to the results 

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

        df['CC'] = df['New_in_species'].apply(is_empty)

        df = pd.merge(df3[['filename']],
                df,
                how= 'left', on='filename')
        df = df.fillna(0) #assumign species not present in LotusDB the number of reported compounds is set to 0
        df.to_csv('../data_out/CC_results.tsv', sep='\t')
        return df
        
    else:
        print ('Similarity Class component not calculated')