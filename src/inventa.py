
import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile



from sklearn.metrics import pairwise_distances
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest
from sklearn import preprocessing
from skbio.stats.ordination import pcoa
from skbio import OrdinationResults

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
    return df

def drop_samples_based_on_string(df,list_of_strings_for_QC_Blank_filter,column):
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

    return df

def reduce_df(df, metadata_df, column):
    """ Reduce the full df to minimal info

    Args:
        df = full_df object (pandas table)

    Returns:
        reduced_df
    """
    reduced_df = df
    reduced_df.set_index(column, inplace=True)
    reduced_df = reduced_df.iloc[:,len(metadata_df.columns)-1:]
    return reduced_df


def top_ions(df1, df2):
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

    #computes the top species/part for each feature 
    df3 = df1.copy().transpose()
    df3 = df3.astype(float)
    df3 = df3.apply(lambda s: s.abs().nlargest(1).index.tolist(), axis=1)
    df3 = df3.to_frame()
    df3[['ATTRIBUTE_Sppart']] = pd.DataFrame(df3[0].values.tolist(),index= df3.index)
    df3 = df3.drop([0], axis=1)
    df3.reset_index(inplace=True)
    df3.rename(columns={'index': 'row ID'}, inplace=True)
    df3['row ID'] = df3['row ID'].astype(int)
    
    #merge all the data 
    df = pd.merge(left=df3,right=dfA, how='left', on='row ID')
    df = pd.merge(left=df2,right=df, how='left',on='row ID')
    return df


def annotations(df1, df2, df3,
                only_ms2_annotations, sirius_annotations, isbd_annotations,
                min_score_final, min_ConfidenceScore, min_ZodiacScore):

    """ 
        function to check the presence of annotations by feature in the combined information form gnps &/ in silico 
        
        Args:
        df1 = clusterinfosummary # mandatory 
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

    if (isbd_annotations == True and sirius_annotations == True):
        #GNPS + ISDB + SIRIUS
        
        #work on gnps annotations

        #find null values (non annotated)
        df1['Annotated'] = pd.isnull(df1['GNPS_INCHI_MF'])
            #lets replace the booleans 
        bD = {True: '0', False: '1'}
        df1['Annotated_GNPS'] = df1['Annotated'].replace(bD)

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
        df = pd.merge(left=df1[['cluster index', 'componentindex', 'Annotated_GNPS']], right=df2[['cluster index','Annotated_ISDB']], 
                        how='left', on= 'cluster index')
        df = pd.merge(left=df, right=df3[['cluster index','Annotated_Sirius']], 
                        how='left',on= 'cluster index')

        def annotations_conditions(df):
            """ function to classify the annotations results 
                Args:
                df = treated and combinend table with the gnps and insilico results
                Returns:
                None
            """
            if (df['Annotated_GNPS'] == '1') | (df['Annotated_ISDB'] == '1') | (df['Annotated_Sirius'] == '1'):
                return 1
            else: 
                return 0

        df['annotation'] = df.apply(annotations_conditions, axis=1)

    elif(isbd_annotations == True and sirius_annotations == False):
        #GNPS + ISDB
        #work on gnps annotations

        #find null values (non annotated)
        df1['Annotated'] = pd.isnull(df1['GNPS_INCHI_MF'])
            #lets replace the booleans 
        bD = {True: '0', False: '1'}
        df1['Annotated_GNPS'] = df1['Annotated'].replace(bD)

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
            
        df = pd.merge(left=df1[['cluster index', 'componentindex', 'Annotated_GNPS']], right=df2[['cluster index','Annotated_ISDB']], 
                        how='left', on= 'cluster index')
        def annotations_conditions(df):
            """ function to classify the annotations results 
                Args:
                df = treated and combinend table with the gnps and insilico results
                Returns:
                None
            """
            if (df['Annotated_GNPS'] == '1') | (df['Annotated_ISDB'] == '1'):
                return 1
            else: 
                return 0

        df['annotation'] = df.apply(annotations_conditions, axis=1) 

    elif(isbd_annotations == False and sirius_annotations == True):
        #GNPS + SIRIUS
        #work on gnps annotations

        #find null values (non annotated)
        df1['Annotated'] = pd.isnull(df1['GNPS_INCHI_MF'])
            #lets replace the booleans 
        bD = {True: '0', False: '1'}
        df1['Annotated_GNPS'] = df1['Annotated'].replace(bD)

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

        df = pd.merge(left=df1[['cluster index', 'componentindex', 'Annotated_GNPS']], right=df3[['cluster index','Annotated_Sirius']], 
                        how='left',on= 'cluster index')

        def annotations_conditions(df):
            """ function to classify the annotations results 
                Args:
                df = treated and combinend table with the gnps and insilico results
                Returns:
                None
            """
            if (df['Annotated_GNPS'] == '1') | (df['Annotated_Sirius'] == '1'):
                return 1
            else: 
                return 0

        df['annotation'] = df.apply(annotations_conditions, axis=1)

    elif(isbd_annotations == False and sirius_annotations == False):
        #ONLY GNPS
        #work on gnps annotations
        
        #find null values (non annotated)
        df1['Annotated'] = pd.isnull(df1['GNPS_INCHI_MF'])
        #lets replace the booleans 
        bD = {True: '0', False: '1'}
        df1['Annotated_GNPS'] = df1['Annotated'].replace(bD)
        #reduced
        df = df1[['cluster index', 'componentindex', 'Annotated_GNPS']]
        df = df.fillna({'Annotated_GNPS':0})

        def annotations_gnps(df):
            """ function to classify the annotations results 
            Args:
            df = treated and combinend table with the gnps and insilico results
            Returns:
            None
            """
            if (df['Annotated_GNPS'] == '1'):
                return 1
            else: 
                return 0

        df['annotation'] = df.apply(annotations_gnps, axis=1)
    else:
        print("INVALID INPUT")

    return df

def feature_component(df1,df2,df3, only_feature_specificity, min_specificity, annotation_preference):
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
    df = pd.merge(df3[['filename', 'ATTRIBUTE_Species', 'ATTRIBUTE_Sppart']], df, how='left', on='filename')
    df = df.sort_values(by=['FC'], ascending=False)
    return df
#literature component
 
def literature_component(df, LC_component, min_comp_reported, max_comp_reported):
    """ function to compute the literature component based on the metadata and combinend information of the Dictionary of natural products and the Lotus DB, 
    Args:
        df2 = metadata_df

    Returns:
        None
    """
    if LC_component == False:
        print('Literature component not calculated')
    else:
        LotusDB = pd.read_csv('../data_loc/LotusDB_inhouse_metadata.csv', 
                       sep=',').dropna()

        #create a set of species from the metadata table
        set_sp = set(df['ATTRIBUTE_Species'].dropna()) #dropna is used to erase all the QC, blanks, etc not having a species associated

        #reduce LotusDB a sp in set 
        LotusDB= LotusDB[LotusDB['organism_taxonomy_09species'].isin(set_sp)]
        LotusDB = LotusDB[['organism_name',
            'organism_taxonomy_06family', 'organism_taxonomy_08genus',
            'organism_taxonomy_09species', 'Reported_comp_Family',
            'Reported_comp_Genus', 'Reported_comp_Species']].drop_duplicates()
        #LotusDB.head()
        
        df = pd.merge(df[['filename', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species']],
                LotusDB[['organism_taxonomy_09species', 'Reported_comp_Family','Reported_comp_Genus', 'Reported_comp_Species']],
                how= 'left', left_on='ATTRIBUTE_Species', right_on='organism_taxonomy_09species').drop_duplicates(subset=['filename'])
        df.drop('organism_taxonomy_09species', axis=1, inplace=True)
        df = df.fillna(0) #assumign species not present in LotusDB the number of reported compounds is set to 0
        df['Reported_comp_Species'] = df['Reported_comp_Species'].astype(int) 

        
        def literature_report(df):
            """ function to give a weigth to the counts of the reported compouds according to the used
            Args:
                df = Literature_component output
            Returns:
                None
            """
            if (df['Reported_comp_Species'] <= min_comp_reported):
                return 1
            elif (df['Reported_comp_Species'] <= max_comp_reported & df['Reported_comp_Species'] >= min_comp_reported): 
                return 0.5
            else:
                return 0

        df['LC'] = df.apply(literature_report, axis=1)
    return df

#similarity component: 

def similarity_component(df, SC_component):
    """ function to compute the similarity component based on the MEMO matrix and machine learning unsupervised clustering methods 
    Args:
        df = memo matrix

    Returns:
        None
    """
    if SC_component == False:
        print('Similarity component not calculated')
    else:
        df1 = df.copy()
        df1.set_index('filename', inplace=True)
        df2 = df.copy()
        
        #specify the parameters of the individual classification algorithms
        clf = IsolationForest(n_estimators=100, 
                    max_samples='auto', 
                    contamination=0.15,
                    max_features=1.0, 
                    bootstrap=False, 
                    n_jobs=None, 
                    random_state=None)
        clf.fit(df1)
        pred = clf.predict(df1)
        df1['anomaly_IF'] = pred
        outliers = df1.loc[df1['anomaly_IF']==-1]
        outlier_index = list(outliers.index)

        lof = LocalOutlierFactor(n_neighbors=10, 
                            algorithm='auto',
                            leaf_size=30,
                            metric='braycurtis', 
                            contamination= 0.15,
                            novelty=False, 
                            n_jobs=None)#-1)
        df1['anomaly_LOF'] = lof.fit_predict(df1)
        outliers = df1.loc[df1['anomaly_LOF']==-1]
        outlier_index = list(outliers.index)

        ocsvm = OneClassSVM(kernel='rbf', 
                            degree=3, 
                            gamma='scale', 
                            tol= 1e-3, 
                            max_iter=-1, 
                            nu=0.01)
        df1['anomaly_OCSVM'] = ocsvm.fit_predict(df1)
        outliers = df1.loc[df1['anomaly_OCSVM']==-1]
        outlier_index = list(outliers.index)

        #recover and print the results
        df1.reset_index(inplace=True)
        df = pd.merge(df1,df2, how='left', left_on='filename', right_on='filename')
        df = df[['filename', 'anomaly_IF', 'anomaly_LOF', 'anomaly_OCSVM']]

        def similarity_conditions(df):
            if (df['anomaly_IF'] == -1) | (df['anomaly_LOF'] == -1) | (df['anomaly_OCSVM'] == -1):
                return 1
            else: 
                return 0 

    df['SC'] = df.apply(similarity_conditions, axis=1)
    return df

#Class component:

def sirius_classes(df1,df2,df3): 
    """ function to find the chemical classes proposed by sirius and assign them to a specific species based on the top specificity of the feature
    Args:
        df1 = specificity_df
        df2 = metadata_df
        df3 = output from SIRIUS + Canopus 

    Returns:
        None
    """
    # merge with top filename with iones 
    df3['shared name'] = df3['name'].str.split('_').str[-1].astype(int)
        #the specificity_df is used to assign the main biological source to each feature. 

    df3 = pd.merge(left=df1[['row ID', 'filename', 'ATTRIBUTE_Sppart']], right=df3[['shared name', 'classe']], how='left', left_on='row ID', right_on='shared name').dropna()
    df3.drop('shared name', axis=1, inplace=True)
    df4 = df3[['filename', 'classe']].groupby('filename').agg(set)
    df = pd.merge(left=df2[['filename', 'ATTRIBUTE_Species']], right=df4, how='left', left_on='filename', right_on='filename').dropna()
    df.drop('ATTRIBUTE_Species', axis=1, inplace=True)
    return df

def search_reported_class(df):
    """ function to search the reported chemical classes in each species of the set 
    Args:
        df = metadata_df
        Returns:
        None
    """
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
    df = pd.merge(df[['filename', 'ATTRIBUTE_Species', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Family', 'ATTRIBUTE_Family', 'ATTRIBUTE_Sppart']],df4,left_on= 'ATTRIBUTE_Species', right_on='organism_taxonomy_09species', how='left')
    df = pd.merge(df,df5,left_on= 'ATTRIBUTE_Genus', right_on='organism_taxonomy_08genus', how='left') 
    return df

def class_component(df1, df2, df3, CC_component):
    """ function to compute the class component based on the possible presence of new chemical classes 
    Args:
        df1 = reported_classes_df 
        df2 = sirius_classes_df
        df3 = metadata_df
        Returns:
        None
    """
    if CC_component == False:
        print ('Similarit Class component not calculated')
    else:
        #merge the both tables
        df = pd.merge(df1,df2,on='filename', how='left').dropna()

        #get the difference between sets 

        df['New_in_species'] = df["classe"] - df["Chemical_class_reported_in_species"]  #check if the chemical classes from Sirius are reported in the species
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
    return df

def priority_rank(df1, df2, df3, df4, LC_component, SC_component, CC_component, w1, w2, w3, w4):
    df = df1
    if LC_component == True: 
        df =pd.merge(
                    left=df,
                    right=df2[['filename', 'Reported_comp_Species', 'Reported_comp_Genus', 'LC', 'ATTRIBUTE_Family']], 
                    how='left', 
                    left_on='filename', 
                    right_on='filename')
    else:
        df

    if SC_component == True:
        df =pd.merge(
                    left=df,
                    right=df3[['filename', 'SC']], 
                    how='left', 
                    left_on='filename', 
                    right_on='filename')
    else:
        df

    if CC_component == True: 
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
    
    return df

def process_gnps_results(gnps_folder_path):
    """ function to compute the class component based on the possible presence of new chemical classes 
    Args:
        gnps_folder_path
        Returns: pandas table (deactivated here) and path
    """

    try:
        path = [x for x in os.listdir(gnps_folder_path+'/result_specnets_DB')][0]
        df_annotations = pd.read_csv(gnps_folder_path+'/result_specnets_DB/'+path, sep='\t')
        print('==================')
        print('Classical MN job detected')
        print('==================')
        print('   Number of spectral library annotations in the job = '+str(df_annotations.shape[0]))
        
        path_networkinfo = [x for x in os.listdir(gnps_folder_path+'/clusterinfosummarygroup_attributes_withIDs_withcomponentID')][0]
        #df_network = pd.read_csv(gnps_folder_path+'/clusterinfosummarygroup_attributes_withIDs_withcomponentID/'+path_networkinfo, sep='\t')
        print('==================')
        print('   Number of network nodes in the job = '+str(df_network.shape[0]))
    
    
    except: 
        print('==================')
        try: 
            path = [x for x in os.listdir(gnps_folder_path+'/DB_result')][0]
            df_annotations = pd.read_csv(gnps_folder_path+'/DB_result/'+path, sep='\t')
            print('==================')
            print('FBMN job detected')
            print('==================')
            print('   Number of spectral library annotations in job = '+str(df_annotations.shape[0]))
            
            path_networkinfo = [x for x in os.listdir(gnps_folder_path+'/clusterinfo_summary')][0]
            clusterinfosummary = gnps_folder_path+'/clusterinfo_summary/'+path_networkinfo
            df_network = pd.read_csv(clusterinfosummary, sep='\t')
            print('==================')
            print('   Number of network nodes in the job = '+str(df_network.shape[0]))
            return clusterinfosummary

        except:
            path = [x for x in os.listdir(gnps_folder_path+'/DB_result')][0]
            df_annotations = pd.read_csv(gnps_folder_path+'/DB_result/'+path, sep='\t')
            print('==================')
            print('FBMN job detected')
            print('==================')
            print('   Number of spectral library annotations in job = '+str(df_annotations.shape[0]))
            
            path_networkinfo = [x for x in os.listdir(gnps_folder_path+'/clusterinfosummarygroup_attributes_withIDs_withcomponentID')][0]
            clusterinfosummary = gnps_folder_path+'/clusterinfosummarygroup_attributes_withIDs_withcomponentID/'+path_networkinfo
            df_network = pd.read_csv(clusterinfosummary, sep='\t')
            print('==================')
            print('   Number of network nodes in the job = '+str(df_network.shape[0]))
            return clusterinfosummary

def get_gnp_db_results(gnps_folder_path):
    """ function to compute the class component based on the possible presence of new chemical classes 
    Args:
        gnps_folder_path
        Returns: pandas table (deactivated here) and path
    """

    try:
        path = [x for x in os.listdir(gnps_folder_path+'/result_specnets_DB')][0]
        df_annotations = pd.read_csv(gnps_folder_path+'/result_specnets_DB/'+path, sep='\t')
        path_networkinfo = [x for x in os.listdir(gnps_folder_path+'/clusterinfosummarygroup_attributes_withIDs_withcomponentID')][0]
     
    except: 
        try: 
            path = [x for x in os.listdir(gnps_folder_path+'/DB_result')][0]
            df_annotations = pd.read_csv(gnps_folder_path+'/DB_result/'+path, sep='\t')
            
            path_networkinfo = [x for x in os.listdir(gnps_folder_path+'/DB_result')][0]
            db_results = gnps_folder_path+'/DB_result/'+path_networkinfo
            df_network = pd.read_csv(db_results, sep='\t')
            return db_results

        except:
            path = [x for x in os.listdir(gnps_folder_path+'/DB_result')][0]
            df_annotations = pd.read_csv(gnps_folder_path+'/DB_result/'+path, sep='\t')
            
            path_networkinfo = [x for x in os.listdir(gnps_folder_path+'/clusterinfosummarygroup_attributes_withIDs_withcomponentID')][0]
            clusterinfosummary = gnps_folder_path+'/clusterinfosummarygroup_attributes_withIDs_withcomponentID/'+path_networkinfo
            df_network = pd.read_csv(clusterinfosummary, sep='\t')
            return clusterinfosummary    