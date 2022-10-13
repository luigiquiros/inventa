
import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib
from tqdm import tqdm


#general treatment 


def features_filter(df, min_threshold):
        
    df[df<min_threshold] = 0 #change all the values lower than x for 0 in the dataframe
    #once the data was filtered, the table is normalized sample-wise
    df = df.apply(lambda x: x/x.max(), axis=0)
    df.to_csv('../data_out/filtered_quant_df.tsv', sep='\t')
    return df

def quantile_filter(df, quantile_threshold):
    
    df = df.replace(0, np.nan)
    df = df.mask(df < df.quantile(quantile_threshold))
    df = df.fillna(0)
    df = df.apply(lambda x: x/x.max(), axis=0)
    return df


def full_data(df1, df2, filename_header):
    """ merge and format the metadata + quantitative information 

    Args:
        df1 = metadata table
        df2 = quantitative.csv file, output from MZmine

    Returns:
        None
    """
    df2 = df2.transpose()
    df2.index.name = filename_header
    df2.reset_index(inplace=True)
    df2.set_index(filename_header, inplace=True)
    df = pd.merge(df1, df2, how='outer', on=filename_header)
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

def drop_samples_based_on_string_ind(repository_path, ionization_mode, filename_header, sampletype_header, metric_df, metadata_df, list_of_strings_for_QC_Blank_filter, column):
    """ drop samples based on string 
    Args:
        pd dataframe
        list of string

    Returns:
        pd dataframe
    """
    #metric_df = pd.read_csv(MEMO_path, sep='\t')
    #metadata_df = pd.read_csv(metadata_path, sep='\t')

    df= pd.merge(metric_df, metadata_df[[filename_header, sampletype_header]], how='left', on=filename_header)
    print(df.shape)

    for string in list_of_strings_for_QC_Blank_filter:
        df = df[~df[column].str.contains(string, na=False)]
        df = df.dropna(how = 'any', subset=[column])
    print(df.shape)
    
    df.drop(sampletype_header, axis=1, inplace=True)
    #df.drop('Unnamed: 0', axis=1, inplace=True)
    path = os.path.normpath(repository_path)
    pathout = os.path.join(path, 'results/')
    os.makedirs(pathout, exist_ok=True)
    pathout = os.path.join(pathout, 'memo_matrix_filtered' +'_' + ionization_mode + '.tsv')
    df.to_csv(pathout, sep ='\t')

    return df
    
def reduce_df(full_df, metadata_df, col_id_unique):
    """ Reduce the full df to minimal info

    Args:
        df = full_df object (pandas table)

    Returns:
            df
    """
    df= full_df
    df.set_index(col_id_unique, inplace=True)
    df = df.iloc[:,len(metadata_df.columns)-1:]
    df.to_csv('../data_out/reduced_df.tsv', sep='\t')
    return df


def priority_rank(FC, LC, SC, CC, LC_component, SC_component, CC_component, w1, w2, w3, w4, filename_header):

    
    if LC_component == True: 

    #df2 = pd.read_csv('../data_out/LC_results.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df =pd.merge(
                left=FC,
                right=LC[[filename_header, 'LC', 'Reported_comp_Species', 'Reported_comp_Genus', 'Reported_comp_Family']], 
                how='left', 
                on=filename_header)
    else:
        df

    if SC_component == True:
        #df3 = pd.read_csv('../data_out/SC_results.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df =pd.merge(
                    left=df,
                    right=SC[[filename_header, 'SC']], 
                    how='left', 
                    on=filename_header)
    else:
        df

    if CC_component == True:
        #df4 = pd.read_csv('../data_out/CC_results.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df =pd.merge( df, right=CC[[filename_header,'CCs','CCg', 'CC', 'New_CC_in_sp', 'New_CC_in_genus']], 
                        how='left', 
                        on =filename_header)
    else: 
        df

    def priority(df):
        df['PS'] = w1*df['FC']
    
        if LC_component == True: 
            df['PS'] = w1*df['FC'] + w2*df['LC']
        else:
            df

        if SC_component == True:
            df['PS'] = w1*df['FC'] + w2*df['LC'] + w3*df['SC']
        else:
            df

        if CC_component == True: 
            df['PS'] = w1*df['FC'] + w2*df['LC'] + w3*df['SC'] + w4*df['CC']
        else: 
            df
        return df

    df = priority(df)
    df.dropna(inplace=True)
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

def selection_changed_FC(selection):
    return FC.iloc[selection]

def selection_changed(selection):
    return PS.iloc[selection]

#Function to count features different from 0 in each sample 
def feature_count(df, header, filename_header):
    '''count total features more than 0 in each sample
    '''
    df = df[df>0.0].count()
    df = pd.DataFrame(df, columns=[header])
    df.reset_index(inplace=True)
    df.rename(columns={'index': filename_header}, inplace=True)
    return df

def priority_score_ind(repository_path, filename_header, ionization_mode, species_column, genus_column, family_column, sppart_column, w1, w2, w3, w4):
    

    path = os.path.normpath(repository_path)
    samples_dir = [directory for directory in os.listdir(path)]

    for directory in tqdm(samples_dir):
        AC_path = os.path.join(path +'/results/', 'Annotation_component_results' +'_'+ ionization_mode + '.tsv')
        LC_path = os.path.join(path +'/results/', 'Literature_component_results.tsv')
        SC_path = os.path.join(path +'/results/', 'Similarity_component_results' +'_'+ ionization_mode + '.tsv')
        CC_path = os.path.join(path +'/results/', 'Class_component_results' +'_'+ ionization_mode + '.tsv')
        metadata_path = os.path.join(path +'/results/', 'Metadata_combined.tsv')

        try:
            df = pd.read_csv(AC_path, sep='\t')
            LC = pd.read_csv(LC_path, sep='\t')
            SC = pd.read_csv(SC_path, sep='\t')
            CC = pd.read_csv(CC_path, sep='\t')
            metadata_df= pd.read_csv(metadata_path, sep='\t')
            
        except FileNotFoundError:
            continue
        except NotADirectoryError:
            continue


    if os.path.isfile(AC_path):
    
            df = pd.read_csv(AC_path, sep='\t')
            df.drop('Unnamed: 0', axis=1, inplace=True)
            df['PS'] = w1*df['AC']

    else:
            metadata_df = pd.read_csv(metadata_path, usecols= [filename_header, family_column, genus_column, species_column, sppart_column], sep='\t')
            df = metadata_df.copy()

    if os.path.isfile(LC_path):
            
            LC = pd.read_csv(LC_path, sep='\t')
            LC.drop('Unnamed: 0', axis=1, inplace=True)
            df =pd.merge(
                left=df,
                right=LC[[filename_header, 'LC', 'Reported_comp_Species', 'Reported_comp_Genus', 'Reported_comp_Family']], 
                how='left', 
                on=filename_header)
            df['PS'] = w1*df['AC'] + w2*df['LC']
    else:
        df

    if os.path.isfile(CC_path):
        CC = pd.read_csv(CC_path, sep='\t')
        CC.drop('Unnamed: 0', axis=1, inplace=True)
        df =pd.merge(
                        left=df,
                        right=CC[[filename_header,'CCs','CCg', 'CC', 'New_CC_in_sp', 'New_CC_in_genus']], 
                        how='left', 
                        on =filename_header)
        df['PS'] = w1*df['AC'] + w2*df['LC'] + w3*df['CC']
    else:
        df
    
    if os.path.isfile(SC_path):
        SC = pd.read_csv(SC_path, sep='\t')
        SC.drop('Unnamed: 0', axis=1, inplace=True)
        df =pd.merge(
                    left=df,
                    right=SC[[filename_header, 'SC']], 
                    how='left', 
                    on=filename_header)
        df['PS'] = w1*df['AC'] + w2*df['LC'] + w3*df['CC'] + w4*df['SC']
    else:
        df
    df = df[['sample_id', 'organism_species', 'organism_genus', 'organism_family',
       'organism_sppart', 'initial_features', 'features_after_filtering',
       'Annot_features_after_filtering', 'AC', 'LC',
       'Reported_comp_Species', 'Reported_comp_Genus', 'Reported_comp_Family',
       'CCs', 'CCg', 'CC', 'New_CC_in_sp', 'New_CC_in_genus', 'SC','PS']]
    path = os.path.normpath(repository_path)
    pathout = os.path.join(path, 'results/')
    os.makedirs(pathout, exist_ok=True)
    pathout = os.path.join(pathout, 'Priority_score_results' +'_'+ ionization_mode + '.tsv')
    df.to_csv(pathout, sep ='\t')
    return df

def selection_changed_AC(selection):
    return AC.iloc[selection]