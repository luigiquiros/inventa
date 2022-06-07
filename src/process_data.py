
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


def quant_table(df, filter = True, min_threshold = 0.5):
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
    # vertical normalization by sample
    df = df.transpose()
    df = df.div(df.sum(axis=1), axis=0)*100
    df = df.transpose()
    df.to_csv('../data_out/quant_df.tsv', sep='\t')
    return df

def features_filter(df, min_threshold):
        
    df[df<min_threshold] = 0 #change all the values lower than x for 0 in the dataframe
    #once the data was filtered, the table is normalized sample-wise
    df = df.transpose()
    df = df.div(df.sum(axis=1), axis=0)*100
    df = df.transpose()
    df.to_csv('../data_out/filtered_quant_df.tsv', sep='\t')
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
    df= pd.read_csv('../data_out/FC_results.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
    
    if LC_component == True: 
        df2 = pd.read_csv('../data_out/LC_results.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df =pd.merge(
                    left=df,
                    right=df2[['filename', 'LC', 'Reported_comp_Species', 'Reported_comp_Genus', 'Reported_comp_Family']], 
                    how='left', 
                    left_on='filename', 
                    right_on='filename')
    else:
        df

    if SC_component == True:
        df3 = pd.read_csv('../data_out/SC_results.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df =pd.merge(
                    left=df,
                    right=df3[['filename', 'SC']], 
                    how='left', 
                    left_on='filename', 
                    right_on='filename')
    else:
        df

    if CC_component == True:
        df4 = pd.read_csv('../data_out/CC_results.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df =pd.merge(
                        left=df,
                        right=df4[['filename','CC', 'New_CC_in_sp', 'New_CC_in_genus']], 
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


def Cyt_format(col_id_unique): 

    #load red dataframe 
    df = pd.read_csv('../data_out/reduced_df.tsv', sep='\t')#.drop(['Unnamed: 0'],axis=1)
    df.set_index(col_id_unique, inplace=True)
    df = df.transpose()

    #normalize values row-wise 
    norm = df.copy()
    norm['Total'] = norm.sum(axis=1)
    norm = norm.drop(norm[norm.Total == 0].index)
    norm = norm.div(norm.Total, axis=0)
    norm = norm*100
    norm.drop('Total', axis=1, inplace=True)
    norm.head()
    df = norm.transpose()

    #load the final PR values 
    PR = pd.read_csv('../data_out/Priority_rank_results.tsv', sep='\t', 
                        usecols =[col_id_unique, 'PR'])#.drop(['Unnamed: 0'],axis=1)

    #merge both df 
    df2 = pd.merge(df, PR, how ='left', on = 'ATTRIBUTE_Sppart')
    df2.set_index(col_id_unique, inplace=True)
    df3 = df2.multiply(df2.PR, axis=0)
    df3.drop('PR', axis=1, inplace=True)
    df3.loc['Score_Total']= df3.sum()

    #recover the usuful info 
    df = df3.transpose()
    df.reset_index(inplace=True)
    df = df[['index','Score_Total']]
    df = df.astype(int)
    df.to_csv('../data_out/PR_cyto_visualization.tsv', sep='\t')
    return df

def quant_plot(df):
    """ Cleans up the quantitative table to specific format

    Args:
        df = quantitative.csv file, output from MZmine

    Returns:
        None
    """
    df.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
    df.rename(columns = lambda x: x.replace('row retention time', 'retention time (min)'),inplace=True)
    df.drop(list(df.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
    #df.drop('row m/z', axis=1, inplace=True)
    #df.drop('row retention time', axis=1, inplace=True)
    #df.to_csv('../data_out/quant_df.tsv', sep='\t')
    return df

#Function to count features different from 0 in each sample 
def feature_count(df, header, filename_header):
    '''count total features more than 0 in each sample
    '''
    df = df[df>0.0].count()
    df = pd.DataFrame(df, columns=[header])
    df.reset_index(inplace=True)
    df.rename(columns={'index': filename_header}, inplace=True)
    return df