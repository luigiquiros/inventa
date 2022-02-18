
import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib

#general treatment 

def quant_table(df):
    """ Cleans up the quantitative table to specific format

    Args:
        df = quantitative.csv file, output from MZmine

    Returns:
        None
    """
    df.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
    df.drop(list(df.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
    df.drop('row m/z', axis=1, inplace=True)
    df.drop('row retention time', axis=1, inplace=True)
    df.to_csv('../data_out/quant_df.tsv', sep='\t')
    return df

def full_data(df1, df2):
    """ merge and format the metadata + quantitative information 

    Args:
        df1 = metadata table
        df2 = quantitative.csv file, output from MZmine

    Returns:
        None
    """
    df2 = df2.transpose()
    df2.index.name = 'filename'
    df2.reset_index(inplace=True)
    df2.set_index('filename', inplace=True)
    df = pd.merge(df1, df2, how='outer', on='filename')
    df.to_csv('../data_out/full_metadata.tsv', sep='\t')
    return df

def drop_samples_based_on_string(df,filename,list_of_strings_for_QC_Blank_filter,column):
    """ drop samples based on string 

    Args:
        pd dataframe
        list of string

    Returns:
        pd dataframe
    """
    print(df.shape)
    for string in list_of_strings_for_QC_Blank_filter:
        df = df[~df[column].str.contains(string, na=False)]
        df = df.dropna(how = 'any', subset=[column])
    print(df.shape)
    save_path = '../data_out/'
    completeName = os.path.join(save_path, filename+".tsv")
    df.to_csv(completeName, sep='\t')
    return df

def reduce_df(col_id_unique):
    """ Reduce the full df to minimal info

    Args:
        df = full_df object (pandas table)

    Returns:
            df
    """
    df= pd.read_csv('../data_out/full_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
    metadata_df= pd.read_csv('../data_out/metadata_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
    df.set_index(col_id_unique, inplace=True)
    df = df.iloc[:,len(metadata_df.columns)-1:]
    df.to_csv('../data_out/reduced_df.tsv', sep='\t')
    return df


def priority_rank(LC_component, SC_component, CC_component, w1, w2, w3, w4):
    df= pd.read_csv('../data_out/FC_results.tsv', sep='\t')
    
    if LC_component == True: 
        df2 = pd.read_csv('../data_out/LC_results.tsv', sep='\t')
        df =pd.merge(
                    left=df,
                    right=df2[['filename', 'Reported_comp_Species', 'Reported_comp_Genus', 'LC', 'ATTRIBUTE_Family']], 
                    how='left', 
                    left_on='filename', 
                    right_on='filename')
    else:
        df

    if SC_component == True:
        df3 = pd.read_csv('../data_out/SC_results.tsv', sep='\t')
        df =pd.merge(
                    left=df,
                    right=df3[['filename', 'SC']], 
                    how='left', 
                    left_on='filename', 
                    right_on='filename')
    else:
        df

    if CC_component == True:
        df2 = pd.read_csv('../data_out/CC_results.tsv', sep='\t') 
        df =pd.merge(
                        left=df,
                        right=df4[['filename', 'New_in_species', 'New_in_genus', 'CC']], 
                        how='left', 
                        left_on='filename', 
                        right_on='filename')
    else: 
        df

    def priority(df):
        df['PR'] = w1*df['FC']
    
        if LC_component == True: 
            df['PR'] = w1*df['FC'] + w2*df['LC']
        else:
            df

        if SC_component == True:
            df['PR'] = w1*df['FC'] + w2*df['LC'] + w3*df['SC']
        else:
            df

        if CC_component == True: 
            df['PR'] = w1*df['FC'] + w2*df['LC'] + w3*df['SC'] + w4*df['CC']
        else: 
            df
        return df

    df = priority(df)
    df.to_csv('../data_out/Priority_rank_results.tsv', sep='\t')
    return df