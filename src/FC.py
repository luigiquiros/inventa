import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib

def top_ions(df1, df2, col_id_unique):
    """ function to compute the top species, top filename and top species/plant part for each ion 
    Args:
        df1 = reduced_df, table of with index on sp/part column and features only.
        df2 = quantitative.csv file, output from MZmine
        Returns:
        None
    """
    #computes the % for each feature
    dfA = df1.copy().transpose()
    dfA = dfA.div(dfA.sum(axis=1), axis=0)*100
    dfA.reset_index(inplace=True)
    dfA.rename(columns={'index': 'row ID'}, inplace=True)
    dfA.set_index('row ID', inplace=True)
    dfA = dfA.astype(float)
    dfA['Feature_specificity'] = dfA.apply(lambda s: s.abs().nlargest(1).sum(), axis=1)
    dfA.reset_index(inplace=True)
    #df1 = df1.drop([0], axis=1)
    dfA = dfA[['row ID', 'Feature_specificity']]

    #computes the top filename for each ion 
    #df2 = quant_df
    df2 = df2.div(df2.sum(axis=1), axis=0)*100
    df2 = df2.copy()
    df2 = df2.astype(float)
    df2 = df2.apply(lambda s: s.abs().nlargest(1).index.tolist(), axis=1)
    df2 = df2.to_frame()
    df2['filename'] = pd.DataFrame(df2[0].values.tolist(), index= df2.index)
    df2 = df2.drop([0], axis=1)

    df = pd.merge(left=dfA,right=df2, how='left',on='row ID')

    if col_id_unique != 'filename':
        #computes the top species/part for each feature 
        df3 = df1.copy().transpose()
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


def annotations(df1, df2, df3,
                only_ms2_annotations, sirius_annotations, isbd_annotations,
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
    df =df1
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
                    
        if only_ms2_annotations == True:
            df2.loc[df2['lib_type']== 'MS1_match', 'Annotated_ISDB'] = 0
        else:
            df2
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
 
    return df 

def feature_component(df1,df2,df3, only_feature_specificity, min_specificity, annotation_preference, col_id_unique):
    """ function to calculate the feature specificity and feature component, as default both columns are added. 
    Args:
        df1 = specificity_df, calculated with the top_ions function 
        df2 = annotation_df, calculated with the annotations function
        df3 = metadata_df
    Returns:
        None
    """
    df4 = pd.merge(df1,df2, how='left', left_on='row ID', right_on='cluster index')

    if only_feature_specificity == True:
        #Computation of the general specificity 
        df5 = df4.copy().groupby('filename').apply(lambda x: len(x[(x['Feature_specificity']>= min_specificity)])).sort_values(ascending=False)
        df5 = df5.div(df4.groupby('filename').Feature_specificity.count(), axis=0)
        df = pd.DataFrame(df5)
        df.rename(columns={0: 'Sample_specificity'}, inplace=True)

    else: 
    #Computation of the general specificity 
        df5 = df4.copy().groupby('filename').apply(lambda x: len(x[(x['Feature_specificity']>= min_specificity)])).sort_values(ascending=False)
        df5 = df5.div(df4.groupby('filename').Feature_specificity.count(), axis=0)
        df5 = pd.DataFrame(df5)
        df5.rename(columns={0: 'Sample_specificity'}, inplace=True)

    #Computation of the feature component 
    df6 = df4.copy().groupby('filename').apply(lambda x: len(x[(x['Feature_specificity']>= min_specificity) & (x['annotation']== annotation_preference)])).sort_values(ascending=False)
    df6 = df6.div(df4.groupby('filename').Feature_specificity.count(), axis=0)
    df6 = pd.DataFrame(df6)
    df6.rename(columns={0: 'FC'}, inplace=True)
    df = pd.merge(df5, df6, how='left', on='filename')

    if col_id_unique != 'filename':
        df = pd.merge(df3[['filename', 'ATTRIBUTE_Species', col_id_unique]], df, how='left', on='filename')
    else:
        df
    df = df.sort_values(by=['FC'], ascending=False)
    df.to_csv('../data_out/FC_results.tsv', sep='\t')
    return df