import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib

def top_ions(col_id_unique):
    """ function to compute the top species, top filename and top species/plant part for each ion 
    Args:
        df1 = reduced_df, table of with index on sp/part column and features only.
        df2 = quantitative.csv file, output from MZmine
        Returns:
        None
    """
    #computes the % for each feature
    dfA = pd.read_csv('../data_out/reduced_df.tsv', sep='\t', index_col=[0])
    dfA = dfA.copy().transpose()
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
    df2 = pd.read_csv('../data_out/quant_df.tsv', sep='\t', index_col=[0])
    df2 = df2.div(df2.sum(axis=1), axis=0)
    df2 = df2.copy()
    df2 = df2.astype(float)
    df2 = df2.apply(lambda s: s.abs().nlargest(1).index.tolist(), axis=1)
    df2 = df2.to_frame()
    df2['filename'] = pd.DataFrame(df2[0].values.tolist(), index= df2.index)
    df2 = df2.drop([0], axis=1)

    df = pd.merge(left=dfA,right=df2, how='left',on='row ID')

    if col_id_unique != 'filename':
        #computes the top species/part for each feature 
        df3 = pd.read_csv('../data_out/reduced_df.tsv', sep='\t', index_col=[0])
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
    df.to_csv('../data_out/specificity_df.tsv', sep='\t')
    return df


def annotations(df2, df3,
                sirius_annotations, isbd_annotations,
                min_score_final, min_ConfidenceScore, min_ZodiacScore):

    """ 
        function to check the presence of annotations by feature in the combined information form gnps &/ in silico 
        
        Args:
        df1 = annot_gnps_df # mandatory 
        df2 = tima_results_filename
        df3 = sirius_annotations_filename
        
        only_ms2_annotations =
        sirius_annotations = 
        isbd_annotations = 
        min_score_final =
        min_ConfidenceScore = 
        min_ZodiacScore  =

        Returns:
        None
    """
    #ONLY GNPS        
    #find null values (non annotated)
    df1 = pd.read_csv('../data_out/annot_gnps_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
    df = df1.copy()
    df['Annotated'] = pd.isnull(df['Consol_InChI'])
    #lets replace the booleans 
    bD = {True: '0', False: '1'}
    df['Annotated_GNPS'] = df['Annotated'].replace(bD)
    #reduced
    df = df[['cluster index', 'componentindex', 'Annotated_GNPS']]
    df = df.fillna({'Annotated_GNPS':0})

    if isbd_annotations == True:
        # work on df2 (isdb annotations)

        df2 = pd.merge(left=df1[['cluster index']], 
                        right=df2, 
                        how='left', left_on= 'cluster index', right_on='feature_id')

        #recover one value from multiple options:
        df2['score_final'] = df2['score_final'].str.split('|').str[-1].astype(float)
        df2['lib_type'] = df2['score_initialNormalized'].str.split('|').str[-1].astype(float)
        df2.drop('score_initialNormalized', axis=1, inplace=True)
        df2['molecular_formula'] = df2['molecular_formula'].str.split('|').str[-1].astype(str)
        
        def score_final_isdb(final_score):
            if final_score >= min_score_final:
                annotated=1 #good annotation
            else:
                annotated=0 #'bad annotation'
            return annotated   

        df2['Annotated_ISDB'] = df2.apply(lambda x: score_final_isdb(x['score_final']), axis=1)
        df2.loc[df2['lib_type']== 'MS1_match', 'Annotated_ISDB'] = 0
     
        #merge the information 
        df = pd.merge(left=df, right=df2[['cluster index','Annotated_ISDB']], 
                        how='left', on= 'cluster index')    
    else:
        df

    if sirius_annotations == True:
        # work on df3 (sirius annotations)
        
        #get the feature id 
        df3['shared name'] = df3['id'].str.split('_').str[-1].astype(int)
        df3 = pd.merge(left=df1[['cluster index']], 
                    right=df3[['shared name','ConfidenceScore','ZodiacScore']], 
                    how='left', left_on= 'cluster index', right_on='shared name')

        df3['ConfidenceScore'] = df3['ConfidenceScore'].fillna(0)

        def Sirius_annotation(ConfidenceScore, ZodiacScore):
            if ConfidenceScore >= min_ConfidenceScore and ZodiacScore >= min_ZodiacScore:
                annotated=1 #good annotation
            else:
                annotated=0 #'bad annotation'
            return annotated

        df3['Annotated_Sirius'] = df3.apply(lambda x: Sirius_annotation(x['ConfidenceScore'], x['ZodiacScore']), axis=1)
            #df3.head(2)
        #merge the information 
        df = pd.merge(left=df, right=df3[['cluster index','Annotated_Sirius']], 
                        how='left',on= 'cluster index')
    else:
        df

    def annotations_gnps(df):
            """ function to classify the annotations results 
            Args:
            df = treated and combinend table with the gnps and insilico results
            Returns:
            None
            """
            if isbd_annotations == True and sirius_annotations == True:

                if (df['Annotated_GNPS'] == '1') | (df['Annotated_ISDB'] == '1') | (df['Annotated_Sirius'] == '1'):
                    return 1
                else: 
                    return 0

            elif isbd_annotations == True and sirius_annotations == False: 
                
                if (df['Annotated_GNPS'] == '1') | (df['Annotated_ISDB'] == '1'):
                    return 1
                else: 
                    return 0

            elif isbd_annotations == False and sirius_annotations == True: 

                if (df['Annotated_GNPS'] == '1') | (df['Annotated_Sirius'] == '1'):
                    return 1
                else: 
                    return 0
            else: 
                if (df['Annotated_GNPS'] == '1'):
                    return 1
                else: 
                    return 0
    df['annotation'] = df.apply(annotations_gnps, axis=1)  
    df.to_csv('../data_out/annotations_df.tsv', sep='\t')
    return df 



def feature_component(quant_df, filtered_quant_df, specificity_df, annotation_df, metadata_df, family_column, genus_column, species_column, col_id_unique, min_specificity, annotation_preference, filename_header, annot_sirius_df, sirius_annotations, annot_gnps_df, min_ZodiacScore):
      
    #1) Feature count by sample before and after filtering
    
    def feature_count(df, header, filename_header):
        df =df[df>0.0].count()
        df = pd.DataFrame(df, columns=[header])
        df.reset_index(inplace=True)
        df.rename(columns={'index': filename_header}, inplace=True)
        return df
    
    initial_features_count = feature_count(quant_df, header ='initial_F', filename_header = filename_header)

    filtered_features_count = feature_count(filtered_quant_df,  header ='filtered_F', filename_header = filename_header)


    #2) combine information from specificity and annotation status for each feature
    df1 = specificity_df.copy()
    df2 = annotation_df.copy()
    df4 = pd.merge(df1,df2, how='left', left_on='row ID', right_on='cluster index')
    
    #3) calculate the total number of features > min_specificity

    specific_features_count = df4.copy().groupby('filename').apply(lambda x: len(x[(x['Feature_specificity']>= min_specificity)])).sort_values(ascending=False)
    specific_features_count = pd.DataFrame(specific_features_count, columns=['Total_SF'])
    specific_features_count.reset_index(inplace=True)
    specific_features_count.rename(columns={'index': filename_header}, inplace=True)

    #4) calculate the total number of features > min_specificity & NON annotated

    specific_non_annotated_features_count = df4.copy().groupby('filename').apply(lambda x: len(x[(x['Feature_specificity']>= min_specificity) & (x['annotation']== annotation_preference)])).sort_values(ascending=False)
    specific_non_annotated_features_count = pd.DataFrame(specific_non_annotated_features_count, columns=['Total_SNAF'])
    specific_non_annotated_features_count.reset_index(inplace=True)
    specific_non_annotated_features_count.rename(columns={'index': filename_header}, inplace=True)

    #5) merge the information
     
    df = pd.merge(initial_features_count, filtered_features_count, how='outer', on=filename_header)
    df = pd.merge(df, specific_features_count, how='outer', on=filename_header)
    df = pd.merge(df, specific_non_annotated_features_count, how='outer', on=filename_header)
    #df.head()

    #6) calculate the ratios
    
    df['FS'] = df['Total_SF']/df['filtered_F']
    df['FS'] =df['FS'].round(decimals = 2)
    df['FC'] = df['Total_SNAF']/df['filtered_F']
    df['FC'] = df['FC'].round(decimals = 2)
    df = df.sort_values(by=['FC'], ascending=False)

    if sirius_annotations == True: 
        df1 = annotation_df.copy() #pd.read_csv('../data_out/annot_gnps_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df2 = annot_sirius_df.copy()
        df2['shared name'] = df2['id'].str.split('_').str[-1].astype(int)
        df3 = specificity_df.copy() #pd.read_csv('../data_out/specificity_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df4 = annotation_df.copy() #pd.read_csv('../data_out/annotations_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df5 = pd.merge(left=df1[['cluster index']],right=df2[['shared name','ZodiacScore']], how='left', left_on= 'cluster index', right_on='shared name')
        df5 = pd.merge(df3,df5, how='left', left_on='row ID', right_on='cluster index')
        df5 = pd.merge(df4[['cluster index','annotation']],df5, how='left', on='cluster index')
        df5.drop('shared name', axis=1, inplace=True)
        df5.drop('cluster index', axis=1, inplace=True)
        df5['ZodiacScore'] = df5['ZodiacScore'].fillna(0)
        
        #calculate the total of features with a good quality MF 
        GQMF_features_count = df5.copy().groupby('filename').apply(lambda x: len(x[(x['Feature_specificity']>= min_specificity) & (x['ZodiacScore'] >= min_ZodiacScore) & (x['annotation'] == annotation_preference)]))
        GQMF_features_count = pd.DataFrame(GQMF_features_count, columns=['Total_SNA_GQMFF'])
        GQMF_features_count.reset_index(inplace=True)
        GQMF_features_count.rename(columns={'index': filename_header}, inplace=True)

        # calculate ratio
        #6) calculate the ratios
        df = pd.merge(df, GQMF_features_count, how='outer', on=filename_header)
        df['MF_prediction_ratio'] = df['Total_SNA_GQMFF']/df['filtered_F']
        df['MF_prediction_ratio'] =df['MF_prediction_ratio'].round(decimals = 2)

        df = df[[filename_header,'initial_F', 'filtered_F', 'Total_SF', 'Total_SNAF', 'Total_SNA_GQMFF', 'MF_prediction_ratio', 'FS', 'FC']]
        df = df.sort_values(by=['FC'], ascending=False)
    else:
        df = df.sort_values(by=['FC'], ascending=False)

    if col_id_unique != filename_header:
        df = pd.merge(metadata_df[[filename_header,family_column, genus_column, species_column, col_id_unique]], df, how='left', on=filename_header)
        df = df.sort_values(by=['FC'], ascending=False)
    else:
        df = pd.merge(metadata_df[[filename_header,family_column, genus_column, species_column]], df, how='left', on=filename_header)
        df = df.sort_values(by=['FC'], ascending=False)


    return df