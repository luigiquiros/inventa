
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

def class_component(df3, min_class_confidence, min_recurrence, CC_component):
    """ function to compute the class component based on the possible presence of new chemical classes 
    Args:
        df1 = reported_classes_df 
        df2 = sirius_classes_df
        df3 = metadata_df
        Returns:
        None
    """
    if CC_component == True:
        df1 = pd.read_csv('../data_out/specificity_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df2 = pd.read_csv('../data_out/metadata_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df3 = canopus_npc_df
        df6 = pd.read_csv('../data_out/metadata_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv',sep=',').dropna()
        # merge with top filename with iones 
        #df3['shared name'] = df3['name'].str.split('_').str[-1].astype(int) #use this line if you don't have a column with the name/shared name

        ## Retrieve classes predicted by Sirius         
        #the specificity_df is used to assign the main biological source to each feature. 
        df3.rename(columns={'name': 'shared name'}, inplace=True)
        df3 = pd.merge(left=df1[['row ID', 'filename']], right=df3[['shared name', 'class', 'classProbability']], how='left', left_on='row ID', right_on='shared name').dropna()
        df3.drop('shared name', axis=1, inplace=True)

        #filter based on min_class_probability
        df3.drop(df3[df3.classProbability > min_class_confidence].index, inplace=True)

        #calculare recurrence of each class 
        df4= df3[['filename', 'class']].groupby(['filename','class']).size().reset_index()
        df4.rename(columns={0: 'recurrence'}, inplace=True)

        df4 = df4[df4['recurrence'] >= min_recurrence].groupby('filename').agg(set)
        df4.drop ('recurrence', axis=1, inplace=True)
                
        df5 = pd.merge(left=df2[['filename', 'ATTRIBUTE_Species']], right=df4, how='left', left_on='filename', right_on='filename').dropna()
        #df5.drop('ATTRIBUTE_Species', axis=1, inplace=True)
        #df5.head(2)
        #df5.to_csv('../data_out/sirus_classes_df.tsv', sep='\t')

        ## Retrieve classes reported in Lotus
            
        #create a set of species present in the metatada and reduce the lotus DB to it
        set_sp = set(df5['ATTRIBUTE_Species'].dropna())
        LotusDB= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]

        #retrieve the chemical classes associated to the species and genus
        df7 = LotusDB[['organism_taxonomy_09species', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_09species').agg(set)
        df7.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_species'}, inplace=True)
        df8 = LotusDB[['organism_taxonomy_08genus', 'structure_taxonomy_npclassifier_03class']].groupby('organism_taxonomy_08genus').agg(set)
        df8.rename(columns={'structure_taxonomy_npclassifier_03class': 'Chemical_class_reported_in_genus'}, inplace=True)

        #merge into a single dataframe
        df9 = pd.merge(df6[['filename', 'ATTRIBUTE_Species', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Family', 'ATTRIBUTE_Family']],df7,left_on= 'ATTRIBUTE_Species', right_on='organism_taxonomy_09species', how='left')
        df9 = pd.merge(df9,df8,left_on= 'ATTRIBUTE_Genus', right_on='organism_taxonomy_08genus', how='left')
        #df9.head(5)
        #df9.to_csv('../data_out/reported_classes_df.tsv', sep='\t')

        #obtain the difference between the predicted and reported compounds
        df = pd.merge(df5,df9,on='filename', how='left').dropna()
        #df10.tail(5)
        #df10.to_csv('../data_out/predicted_and_reported_classes_df.tsv', sep='\t')

        df['New_in_species'] = df["class"] - df["Chemical_class_reported_in_species"]  #check if the chemical classes from Sirius are reported in the species
        df['New_in_genus'] = df["New_in_species"] - df["Chemical_class_reported_in_genus"]

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

        df = pd.merge(df2[['filename']], df,how= 'left', on='filename')
        df = df.fillna(0)
        df.to_csv('../data_out/CC_results.tsv', sep='\t')
        return df
    else:
        print ('No search was done because the Class component is not going to be calculated')