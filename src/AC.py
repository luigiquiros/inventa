from re import T
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

def ind_quant_table_full(repository_path, ionization_mode, file_extention, data_process_origin, use_ion_identity, min_score_final, min_ConfidenceScore, min_ZodiacScore):
        
    path = os.path.normpath(repository_path)
    samples_dir = [directory for directory in os.listdir(path)]

    for directory in tqdm(samples_dir):
        quant_path = os.path.join(path, directory, ionization_mode, directory + '_features_quant_' + ionization_mode + '.csv')
        isdb_path = os.path.join(path, directory, ionization_mode+'/isdb/', directory +'_isdb_reweighted_' + ionization_mode+ '.tsv')
        sirius_path = os.path.join(path, directory, ionization_mode, directory + '_WORKSPACE_SIRIUS', 'compound_identifications_adducts.tsv')
        try:
            df =pd.read_csv(quant_path, sep=',')
            dfs = pd.read_csv(sirius_path, sep='\t')
            dfi = pd.read_csv(isdb_path, sep='\t')
        except FileNotFoundError:
            continue
        except NotADirectoryError:
            continue

        df = pd.read_csv(quant_path, sep= ',')
        df.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
        df.rename(columns = lambda x: x.replace(file_extention, ''),inplace=True)
        df.drop(list(df.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
        df.sort_index(axis=1, inplace=True)

        if data_process_origin == 'MZMine3':


            if use_ion_identity == True:

                df.drop(['row ion mobility',
                    'row ion mobility unit', 'row CCS', 'best ion',
                    'correlation group ID', 'auto MS2 verify',
                    'identified by n=', 'partners', 'neutral M mass'], axis=1, inplace=True)

                #complete correlation groups
                df['annotation network number'] = df['annotation network number'].fillna(df['row ID'].apply(str) + 'x')
                df.drop('row ID', axis =1, inplace=True)
                df = df.groupby('annotation network number', dropna=False).max()

            else:
                #prepare quant table acordingly 

                df.drop(['row ion mobility', 'correlation group ID', 'best ion', 'row ion mobility unit', 'row CCS', 
                'annotation network number', 'auto MS2 verify', 'identified by n=', 'partners', 'neutral M mass'], axis=1, inplace=True)
                df.set_index('row ID', inplace=True)

        else:
            df 

        df1 = df[['row ID', 'row m/z', 'row retention time']]
        df.drop(['row m/z', 'row retention time'], axis =1, inplace=True)

        #normalize
        df.set_index('row ID', inplace=True)
        #df = df.apply(lambda x: x/x.sum(), axis=0)
        df =pd.merge(df, df1, how ='left', on='row ID')

        #rename columns
        df.rename(columns={'row m/z': 'm/z', 'row retention time':'retention time (min)'}, inplace=True)
        df['m/z']=df['m/z'].round(decimals = 6)
        df['retention time (min)']=df['retention time (min)'].round(decimals = 2)

        if os.path.isfile(isdb_path):

            dfi = pd.read_csv(isdb_path, sep='\t', usecols =['feature_id', 'libname', 'structure_molecular_formula','final_score'])

            #recover one value from multiple options:
            dfi['final_score'] = dfi['final_score'].astype(str).str.split('|').str[-1].astype(float)
            dfi['libname'] = dfi['libname'].str.split('|').str[-1].astype(str)
            dfi['structure_molecular_formula'] = dfi['structure_molecular_formula'].str.split('|').str[-1].astype(str)
            dfi['final_score'] = dfi['final_score'].astype(int)

            #quality annotations filtering

            def score_final_isdb(final_score):
                if final_score >= min_score_final:
                    annotated=1 #good annotation
                else:
                    annotated=0 #'bad annotation'
                return annotated   

            dfi['Annotated_ISDB'] = dfi.apply(lambda x: score_final_isdb(x['final_score']), axis=1)
            dfi.drop('final_score', axis =1, inplace=True)
            dfi.loc[dfi['libname']== 'MS1_match', 'Annotated_ISDB'] = 0

            dfi.rename(columns={'libname':'libname_ISDB', 'structure_molecular_formula' : 'structure_molecular_formula_ISDB'}, inplace=True)
            df = pd.merge(df, dfi, how='left', left_on='row ID', right_on='feature_id')
            df.drop('feature_id', axis=1, inplace = True)
        else:
            df

        if os.path.isfile(sirius_path):

            dfs = pd.read_csv(sirius_path, sep='\t', 
                        usecols =['id','molecularFormula', 'InChI', 'ConfidenceScore','ZodiacScore', 'adduct', 'name'])
            
            dfs['shared name'] = dfs['id'].str.split('_').str[-1].astype(int)
            dfs['ConfidenceScore'] = dfs['ConfidenceScore'].fillna(0)
            dfs['ZodiacScore'] = dfs['ZodiacScore'].fillna(0)
            dfs.drop('id', axis=1, inplace = True)
            #df.astype('int64')

            def Sirius_annotation(ConfidenceScore, ZodiacScore):
                if ConfidenceScore >= min_ConfidenceScore and ZodiacScore >= min_ZodiacScore:
                    annotated=1 #good annotation
                else:
                    annotated=0 #'bad annotation'
                return annotated

            dfs['Annotated_Sirius'] = dfs.apply(lambda x: Sirius_annotation(x['ConfidenceScore'], x['ZodiacScore']), axis=1)
            dfs.rename(columns={'molecularFormula':'structure_molecular_formula_SIRIUS', 'adduct': 'adduct_SIRIUS', 'name':'Compound_name_SIRIUS', 'InChI':'structure_inchi_SIRIUS'}, inplace=True)
            dfs.drop(['ConfidenceScore','ZodiacScore'], axis=1, inplace = True)
            df = pd.merge(df, dfs, how='left', left_on='row ID', right_on='shared name')
            df.drop('shared name', axis=1, inplace = True)
        else:
            df
            
        def annotation_status(df):
            
            """ function to classify the annotations results 
            Args:
            df = treated and combinend table with the gnps and insilico results
            Returns:
            None
            """
            if (df['Annotated_ISDB'] == 1) or (df['Annotated_Sirius'] == 1):
                return 1
            else: 
                return 0

        df['annotation'] = df.apply(annotation_status, axis=1)
        pathout = os.path.join(path, 'results/')
        os.makedirs(pathout, exist_ok=True)
        pathout = os.path.join(pathout, directory +'_' + ionization_mode + '_quant_annotations.tsv')
        #pathout = os.path.normpath(os.path.join((pathout, directory +'_' + ionization_mode + '_quant_annotations.tsv')))
        df.to_csv(pathout, sep ='\t')
    print(f'Result are in : {pathout}')

def get_sirius_annotations_ind(repository_path, sirius_sample_suffix, min_ConfidenceScore, min_ZodiacScore):
    df = pd.DataFrame()
    if sirius_annotations == True:
        for r, d, f in os.walk(repository_path):
            for file in (f for f in f if f.endswith(sirius_sample_suffix)):
                    
                    complete_file_path =r+'/'+file 
                    df = pd.read_csv(complete_file_path,sep='\t', 
                                usecols =['id','molecularFormula', 'ConfidenceScore','ZodiacScore', 'adduct', 'name'], 
                                low_memory=False)
                    
                    df['shared name'] = df['id'].str.split('_').str[-1].astype(int)
                    df['ConfidenceScore'] = df['ConfidenceScore'].fillna(0)
                    df['ZodiacScore'] = df['ZodiacScore'].fillna(0)
                    df.drop('id', axis=1, inplace = True)
                    #df.astype('int64')

                    def Sirius_annotation(ConfidenceScore, ZodiacScore):
                        if ConfidenceScore >= min_ConfidenceScore and ZodiacScore >= min_ZodiacScore:
                            annotated=1 #good annotation
                        else:
                            annotated=0 #'bad annotation'
                        return annotated

                    df['Annotated_Sirius'] = df.apply(lambda x: Sirius_annotation(x['ConfidenceScore'], x['ZodiacScore']), axis=1)

                    prefix = 'treated_'
                    df.to_csv(r+'/'+prefix+file, sep ='\t')

def get_isdb_annotations_ind(repository_path, isdb_sample_suffix, isdb_annotations, min_score_final):
    df = pd.DataFrame()
    if isdb_annotations == True:
        for r, d, f in os.walk(repository_path):
            for file in (f for f in f if f.endswith(isdb_sample_suffix)):
                    
                    complete_file_path =r+'/'+file 
                    df = pd.read_csv(complete_file_path, sep='\t', usecols =['feature_id', 'libname', 'structure_molecular_formula','structure_inchi','final_score'], 
                                low_memory=False)
                    
                    #recover one value from multiple options:
                    df['final_score'] = df['final_score'].astype(str).str.split('|').str[-1].astype(float)
                    df['libname'] = df['libname'].str.split('|').str[-1].astype(str)
                    df['structure_molecular_formula'] = df['structure_molecular_formula'].str.split('|').str[-1].astype(str)

                    #quality annotations filtering

                    def score_final_isdb(final_score):
                        if final_score >= min_score_final:
                            annotated=1 #good annotation
                        else:
                            annotated=0 #'bad annotation'
                        return annotated   

                    df['Annotation_ISDB'] = df.apply(lambda x: score_final_isdb(x['final_score']), axis=1)
                    df.drop('final_score', axis =1, inplace=True)
                    df.loc[df['libname']== 'MS1_match', 'Annotated_ISDB'] = 0
    
                    prefix = 'treated_'
                    df.to_csv(r+'/'+prefix+file, sep ='\t')
#AC.X. Load quant information / add annotation status / filter 

def ind_quant_table(repository_path, quant_table_suffix, data_process_origin, use_ion_dentity):

    df = pd.DataFrame()
    for r, d, f in os.walk(repository_path):
        for file in (f for f in f if f.endswith(quant_table_suffix)):
                
                complete_file_path =r+'/'+file 
                df = pd.read_csv(complete_file_path)
                
                df.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
                df.rename(columns = lambda x: x.replace(file_extention, ''),inplace=True)
                df.drop(list(df.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
                df.sort_index(axis=1, inplace=True)

                if data_process_origin == 'MZMine3':


                    if use_ion_dentity == True:

                        df.drop(['row ion mobility',
                            'row ion mobility unit', 'row CCS', 'best ion',
                            'correlation group ID', 'auto MS2 verify',
                            'identified by n=', 'partners', 'neutral M mass'], axis=1, inplace=True)

                        #complete correlation groups
                        df['annotation network number'] = df['annotation network number'].fillna(df['row ID'].apply(str) + 'x')
                        df.drop('row ID', axis =1, inplace=True)
                        df = df.groupby('annotation network number', dropna=False).max()

                    else:
                        #prepare quant table acordingly 

                        df.drop(['row ion mobility', 'correlation group ID', 'best ion', 'row ion mobility unit', 'row CCS', 
                        'annotation network number', 'auto MS2 verify', 'identified by n=', 'partners', 'neutral M mass'], axis=1, inplace=True)
                        df.set_index('row ID', inplace=True)

                else:
                    df 

                df1 = df[['row ID', 'row m/z', 'row retention time']]
                df.drop(['row m/z', 'row retention time'], axis =1, inplace=True)

                #normalize
                df.set_index('row ID', inplace=True)
                df = df.apply(lambda x: x/x.sum(), axis=0)

                df =pd.merge(df, df1, how ='left', on='row ID')

                #rename columns
                df.rename(columns={'row m/z': 'm/z', 'row retention time':'retention time (min)'}, inplace=True)
                df['m/z']=df['m/z'].round(decimals = 6)
                df['retention time (min)']=df['retention time (min)'].round(decimals = 2)

                #add ISDB and Sirius annotations


                #filter


                prefix = 'treated_'
                df.to_csv(r+'/'+prefix+file, sep =',')

def annotation_component(repository_path, ionization_mode, file_extention, intensity_filter, quantile_filter, min_threshold, quantile_threshold):
    
    path = os.path.normpath(repository_path)
    samples_dir = [directory for directory in os.listdir(path)]

    files = []
    original_feature_count = []
    feature_count_filtered = []
    annotated_features_count = []

    for directory in tqdm(samples_dir):
        
        quant_annotations_path = os.path.join(path, path +'/results/', directory + '_'+ionization_mode + '_quant_annotations.tsv')
        
        #get the filename column associated to each quant table
        column = os.path.join(path, directory, directory + '.mzXML')
        column = column.rsplit('/',1)[1]

        try:
            df_original = pd.read_csv(quant_annotations_path, sep='\t')
        except FileNotFoundError:
            continue
        except NotADirectoryError:
            continue
        
        #recover filenames 
        files.append(directory)
        #read original filename
        df_original = pd.read_csv(quant_annotations_path, sep='\t')
        #normalize
        df_original[column] = df_original[column]/df_original[column].sum()

        #original number of features
        dfo = df_original[column]
        dfo = dfo[dfo>0.0].count()
        original_feature_count.append(dfo)

        #check and apply filtering steps if applicable

        if intensity_filter == True and quantile_filter == True:
                
            dff = df_original.copy()
            #apply intensity filter
            dff[column].values[dff[column] < min_threshold] = 0 #change all the values lower than x for 0 in the dataframe
            dff[column] = dff[column]/dff[column].sum() #once the data was filtered, the table is normalized sample-wise

            #apply quantile filtering
            dff = dff.replace (0, np.nan)
            dff = dff[dff[column] < dff[column].quantile(quantile_threshold)]#change all the values lower than x quantile for 0 in the dataframe
            dff[column] = dff[column]/dff[column].sum()#once the data was filtered, the table is normalized sample-wis

        elif intensity_filter == True and quantile_filter == False:
            
            dff = df_original.copy()
            #apply intensity filter
            dff[column].values[dff[column] < min_threshold] = 0 #change all the values lower than x for 0 in the dataframe
            dff[column] = dff[column]/dff[column].sum() #once the data was filtered, the table is normalized sample-wise

        elif intensity_filter == False and quantile_filter == True:
            dff = df_original.copy()
            #apply quantile filtering
            dff = dff.replace (0, np.nan)
            dff = dff[dff[column] < dff[column].quantile(quantile_threshold)]#change all the values lower than x quantile for 0 in the dataframe
            dff[column] = dff[column]/dff[column].sum()#once the data was filtered, the table is normalized sample-wis

        else:
            dff = df_original

        #number of features after filtering
        dffc = dff[column]
        dffc = dffc[dffc>0.0].count()
        feature_count_filtered.append(dffc)


        #number of features after filtering annotated
        dfa = dff[[column, 'annotation']]
        dfa = dfa[dfa['annotation'] == 1]
        dfa = dfa[column]
        dfac = dfa[dfa>0.0].count()
        annotated_features_count.append(dfac)

    AC = pd.DataFrame({filename_header: files,'initial_features': original_feature_count, 'features_after_filtering' : feature_count_filtered, 'Annot_features_after_filtering': annotated_features_count })
    AC['AC'] = (AC['features_after_filtering'] - AC['Annot_features_after_filtering'])/AC['features_after_filtering'] # % of features unannotated
    AC['AC'] = AC['AC'].round(decimals = 1)
    return AC