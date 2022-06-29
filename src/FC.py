import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib

def annotations(df2, df3,
                sirius_annotations, isbd_annotations,
                min_score_final, min_ConfidenceScore, min_ZodiacScore,  correlation_groups_df, use_ion_identity):

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
        df3['ZodiacScore'] = df3['ZodiacScore'].fillna(0)
        df3.drop('shared name', axis=1, inplace = True)
        df3.astype('int64')

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
    df.rename(columns={'cluster index': 'row ID'}, inplace=True)
    
    if use_ion_identity == True: 
        df = pd.merge(correlation_groups_df[['row ID', 'annotation network number']], df[['row ID', 'annotation']], how ='left', on='row ID')
        df.drop('row ID', axis = 1, inplace=True)
        df = df.groupby('annotation network number', as_index=False).agg(max)
        #df.reset_index(inplace=True)
        #df.rename(columns={'index': filename_header}, inplace=True)
    else: 
        df
    #df.to_csv('../data_out/annotations_df.tsv', sep='\t')
    return df 


def feature_component(quant_df, reduced_df, annotation_df, metadata_df, family_column, genus_column, species_column, col_id_unique, min_specificity, annotation_preference, filename_header, annot_sirius_df, sirius_annotations, annot_gnps_df, min_ZodiacScore, multiple_organism_parts, max_parts_per_organism, use_ion_identity):
      
    #1) Feature count by sample before and after filtering
    
    def feature_count(df, header, filename_header):
        df =df[df>0.0].count()
        df = pd.DataFrame(df, columns=[header])
        df.reset_index(inplace=True)
        df.rename(columns={'index': filename_header}, inplace=True)
        return df
    
    #get the number of features > 0 for each sample
    initial_features_count = feature_count(quant_df, header ='initial_F', filename_header = filename_header)
    #get the number of features > min_specificity for each sample
    filtered_features_count = feature_count(reduced_df,  header ='filtered_F', filename_header = filename_header)

    if use_ion_identity == True: 
        row_ID_header = 'annotation network number'
    else: 
        row_ID_header = 'row ID'

    #2) normalize the filtered table and combine information from specificity and annotation status for each feature
    #normalize row-wise the area of features = relative % of each feature in each sample
    
    filtered_quant_df_norm = reduced_df.copy()#.transpose()
    filtered_quant_df_norm = filtered_quant_df_norm.div(filtered_quant_df_norm.sum(axis=1), axis=0).fillna(0)
    filtered_quant_df_norm.reset_index(inplace=True)
    filtered_quant_df_norm.rename(columns={'index': row_ID_header}, inplace=True)
    filtered_quant_df_norm.set_index(row_ID_header, inplace=True)
    filtered_quant_df_norm = filtered_quant_df_norm.astype(float)

    #filtered_quant_df_norm.head(2)
    
    if multiple_organism_parts == True:
        filtered_quant_df_norm['nlargestsum'] = filtered_quant_df_norm.apply(lambda s: s.abs().nlargest(max_parts_per_organism).sum(), axis=1)
        #df.head(2)
        df1 = filtered_quant_df_norm.iloc[:,:-1]
        #compare for greater or equal by division nlargestsum with N and if match replace values
        filtered_quant_df_norm.update(df1.mask(df1.ge(filtered_quant_df_norm['nlargestsum'].div(max_parts_per_organism), axis=0), filtered_quant_df_norm['nlargestsum'], axis=0))
        filtered_quant_df_norm.drop('nlargestsum', axis=1, inplace=True)

        #add the row ID and annotation status of each ion
        filtered_quant_df_norm = pd.merge(annotation_df[[row_ID_header, 'annotation']], filtered_quant_df_norm, how='left', left_on=row_ID_header, right_on=row_ID_header).fillna(0)
        #filtered_quant_df_norm.rename(columns={'cluster index': row_ID_header}, inplace=True)
    else:
        #add the row ID and annotation status of each ion
        filtered_quant_df_norm = pd.merge(annotation_df[[row_ID_header, 'annotation']], filtered_quant_df_norm, how='left', left_on=row_ID_header, right_on=row_ID_header).fillna(0)
        filtered_quant_df_norm.rename(columns={'row ID': row_ID_header}, inplace=True)

    
    #3) calculate the total number of features > min_specificity

        #to get the number of features with area > min_specificity : 
        #select only sample columns
        #check if values in columns >= min_specificity
        #sum all (boolean) values
        #transpose your dataframe

    specific_features_count = filtered_quant_df_norm.copy()
    specific_features_count = specific_features_count.iloc[:,2:].ge(min_specificity).sum().T
    specific_features_count = pd.DataFrame(specific_features_count, columns=['Total_SF'])
    specific_features_count.reset_index(inplace=True)
    specific_features_count.rename(columns={'index': filename_header}, inplace=True)

    #4) calculate the total number of features > min_specificity & NON annotated
        #to get the number of features with area > min_specificity and annotation == annotation preference: 
        #keep only rows where annotation == annotation_preference
        #select only sample columns
        #check if values in columns >= min_specificity
        #sum all (boolean) values
        #transpose your dataframe
    specific_non_annotated_features_count = filtered_quant_df_norm.copy()
    specific_non_annotated_features_count = specific_non_annotated_features_count[specific_non_annotated_features_count['annotation']==annotation_preference].iloc[:,2:].ge(min_specificity).sum().T
    specific_non_annotated_features_count = pd.DataFrame(specific_non_annotated_features_count, columns=['Total_SNAF'])
    specific_non_annotated_features_count.reset_index(inplace=True)
    specific_non_annotated_features_count.rename(columns={'index': filename_header}, inplace=True)
    specific_non_annotated_features_count

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

    if sirius_annotations == True and use_ion_identity == False: 
        df1 = annotation_df.copy() #pd.read_csv('../data_out/annot_gnps_df.tsv', sep='\t').drop(['Unnamed: 0'],axis=1)
        df2 = annot_sirius_df.copy()
        df2['shared name'] = df2['id'].str.split('_').str[-1].astype(int)
        df5 = pd.merge(left=df1[[row_ID_header]],right=df2[['shared name','ZodiacScore']], how='left', left_on= row_ID_header, right_on='shared name')
        df5 = pd.merge( df5, filtered_quant_df_norm, how='left', left_on=row_ID_header, right_on=row_ID_header)
        df5.drop('shared name', axis=1, inplace=True)
        df5.drop(row_ID_header, axis=1, inplace=True)
        df5['ZodiacScore'] = df5['ZodiacScore'].fillna(0)

        #to get the number of features with area > min_specificity and annotation == annotation preference: 
                #keep only rows where annotation == annotation_preference & ['ZodiacScore'] >= min_ZodiacScore
                #select only sample columns
                #check if values in columns >= min_specificity
                #sum all (boolean) values
                #transpose your dataframe

        GQMF_features_count = df5.copy()
        GQMF_features_count = GQMF_features_count[GQMF_features_count['annotation']==annotation_preference]
        GQMF_features_count = GQMF_features_count[GQMF_features_count['ZodiacScore']>=min_ZodiacScore].iloc[:,2:].ge(min_specificity).sum().T
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
    df.to_csv('../data_out/FC_results.tsv', sep='\t')

    return df