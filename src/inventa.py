
#Feature component
FC_component = True  #FC will be calculated
min_specificity = 90  #minimun feature specificity to consider
only_feature_specificity = False  #True if annotations should be ignore and the FC should be calculated based on the features specificity. If False it will compute both The Sample specifity adn the FC
only_gnps_annotations = False  #only the annotations from gpns will be considered 
only_ms2_annotations = False  #False to considere both, MS1 & MS2 annotations, False will only considerer MS2 annotations
annotation_preference= 0  #Only Annotated nodes: '1' 
                          #Only Not annotated: '0'

#Literature component 
LC_component = True  #LC will be calculated
max_comp_reported = 40  #more than this value, the plant is considered no interesting LC =0
min_comp_reported = 10  #less than this value, the plant is consireded very interesintg LC =1
                        #a sample with x between both values gets a LC=0.5

#Class component+
CC_component = True  #CC will be calculated

#Similarity component
SC_component = True  #SC will be calculated

# dependencies 

import pandas as pd
import numpy as np
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
    df.set_index('row ID', inplace=True)
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

def annotations(df1, df2):
    """ function to check the presence of annotations by feature in the combined information form gnps &/ in silico 
    Args:
        df1 = cluster summary results file from GNPS
        df2 = in silico dereplication results file
        Returns:
        None
    """
    if only_gnps_annotations == True:
        #work on gnps annotations
        #df1 = annot_df.copy()
        #find null values (non annotated)
        df1['Annotated'] = pd.isnull(df1['SpectrumID'])
        #lets replace the booleans 
        bD = {True: '0', False: '1'}
        df1['Annotated_GNPS'] = df1['Annotated'].replace(bD)
        #merge the information 
        df = df1[['cluster index', 'componentindex', 'Annotated_GNPS']]
        return df
    else:
        #work on gnps annotations
         #find null values (non annotated)
        df1['Annotated'] = pd.isnull(df1['SpectrumID'])
        #lets replace the booleans 
        bD = {True: '0', False: '1'}
        df1['Annotated_GNPS'] = df1['Annotated'].replace(bD)

        #work on isdb annotations
        if only_ms2_annotations == True:
            df2 = df2[~df2.libname.str.contains('MS1_match', na=False)]
            df2['Annotated'] = pd.isnull(df2['short_inchikey'])
            df2['Annotated_ISDB'] = df2['Annotated'].replace(bD)
        else:
            df2['Annotated'] = pd.isnull(df2['short_inchikey'])
            df2['Annotated_ISDB'] = df2['Annotated'].replace(bD)
        
        #merge the information 
    df = pd.merge(left=df1[['cluster index', 'componentindex', 'Annotated_GNPS']], 
                  right=df2[['feature_id','Annotated_ISDB']], 
                  how='left', left_on= 'cluster index', right_on='feature_id')
    df.drop('feature_id', axis=1, inplace=True)
    df = df.fillna({'Annotated_ISDB':0})
    df['annotation'] = df.apply(annotations_conditions, axis=1)
    return df

def feature_component(df1,df2,df3):
    """ function to calculate the feature specificity and feature component, as default both columns are added. 
    Args:
        df1 = specificity_df, calculated with the top_ions function 
        df2 = annotation_df, calculated with the annotations function
        df3 = metadata_df
    Returns:
        None
    """
    if FC_component == False:
        print('Feature component not calculated')
    else:
        df4 = pd.merge(df1,df2, how='left', left_on='row ID', right_on='cluster index')

        if only_feature_specificity == True:
            #Computation of the general specificity 
            df5 = df4.copy().groupby('filename').apply(lambda x: len(
                x[(x['Feature_specificity']>= min_specificity)])).sort_values(ascending=False)
            df5 = df5.div(df4.groupby('filename').Feature_specificity.count(), axis=0)
            df = pd.DataFrame(df5)
            df.rename(columns={0: 'Sample_specificity'}, inplace=True)

        else: 
            #Computation of the general specificity 
            df5 = df4.copy().groupby('filename').apply(lambda x: len(
                x[(x['Feature_specificity']>= min_specificity)])).sort_values(ascending=False)
            df5 = df5.div(df4.groupby('filename').Feature_specificity.count(), axis=0)
            df5 = pd.DataFrame(df5)
            df5.rename(columns={0: 'Sample_specificity'}, inplace=True)

            #Computation of the feature component 

            df6 = df4.copy().groupby('filename').apply(lambda x: len(
                x[(x['Feature_specificity']>= min_specificity) & (x['annotation']== annotation_preference)])).sort_values(ascending=False)
            df6 = df6.div(df4.groupby('filename').Feature_specificity.count(), axis=0)
            df6 = pd.DataFrame(df6)
            df6.rename(columns={0: 'FC'}, inplace=True)
            df = pd.merge(df5, df6, how='left', on='filename')

        df = pd.merge(df3[['filename', 'ATTRIBUTE_Species', 'ATTRIBUTE_Sppart']], df, how='left', on='filename')
        df = df.sort_values(by=['FC'], ascending=False)
        return df

#literature component
 
def literature_report(y):
    """ function to give a weigth to the counts of the reported compouds according to the used
    Args:
        df1 = Literature_component output
    Returns:
        None
    """
    if (y['Reported_comp_Species'] <= min_comp_reported):
        return 1
    elif (y['Reported_comp_Species'] <= max_comp_reported & y['Reported_comp_Species'] >= min_comp_reported): 
        return 0.5
    else:
        return 0

def literature_component(df):
    """ function to compute the literature component based on the metadata and combinend information of the Dictionary of natural products and the Lotus DB, 
    Args:
        df2 = metadata_df

    Returns:
        None
    """
    if LC_component == False:
        print('Literature component not calculated')
    else:
        LotusDB = pd.read_csv('../data/LotusDB_inhouse_metadata.csv', 
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
        
        df = pd.merge(df2[['filename', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species']],
                LotusDB[['organism_taxonomy_09species', 'Reported_comp_Family','Reported_comp_Genus', 'Reported_comp_Species']],
                how= 'left', left_on='ATTRIBUTE_Species', right_on='organism_taxonomy_09species')

        df = df.fillna(0) #assumign species not present in LotusDB the number of reported compounds is set to 0
        df['Reported_comp_Species'] = df['Reported_comp_Species'].astype(int) 
        df['LC'] = df.apply(literature_report, axis=1)
        return df

#similarity component: 

def similarity_conditions(df):
    if (df['anomaly_IF'] == -1) | (df['anomaly_LOF'] == -1) | (df['anomaly_OCSVM'] == -1):
        return 1
    else: 
        return 0 

def similarity_component(df):
    """ function to compute the similarity component based on the MEMO matrix and machine learning unsupervised clustering methods 
    Args:
        df = meme matrix

    Returns:
        None
    """
    if SC_component == False:
        print('Similarity component not calculated')
    else:
        df2 = df.copy()
        df2.rename(columns = {'Unnamed: 0': 'filename'}, inplace=True)
        columns_to_model=df2.columns[1:39121] #specify the X metrics column names to be modelled (THIS CORRESPOND TO THE SIZE OF THE FEATURES)
        df1 = df2[columns_to_model].astype(np.uint8)
    
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
        df2.reset_index(inplace=True)
        df = pd.merge(df1,df2, how='left', left_on='index', right_on='index')
        df = df[['filename', 'anomaly_IF', 'anomaly_LOF', 'anomaly_OCSVM']]
        df['SC'] = df.apply(similarity_conditions, axis=1)
        return df

#Class component:

def sirius_classes(df1,df2,df3): 
    """ function to find the chemical classes proposed by sirius and assign them to a specific species based on the top specificity of the feature
    Args:
        df1 = specificity_df
        df2 = metadata_df
        df3 = output from SIRIUS and Canopus 

    Returns:
        None
    """
    # merge with top filename with iones 
    df3 = pd.merge(left=df1[['row ID', 'filename', 'ATTRIBUTE_Sppart']], right=df3, how='left', left_on='row ID', right_on='shared name').dropna()
    df3.drop('shared name', axis=1, inplace=True)
    df4 = df3[['filename', 'CAN_classe']].groupby('filename').agg(set)
    df = pd.merge(left=df2[['filename', 'ATTRIBUTE_Species']], right=df4, how='left', left_on='filename', right_on='filename').dropna()
    df.drop('ATTRIBUTE_Species', axis=1, inplace=True)
    return df

def search_reported_class(df):
    """ function to search the reported chemical classes in each species of the set 
    Args:
        df2 = metadata_df
        Returns:
        None
    """
    LotusDB = pd.read_csv('../data/dataLotusDB_inhouse_metadata.csv',
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
    df = pd.merge(df[['filename', 'ATTRIBUTE_Species', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Sppart']],df4,left_on= 'ATTRIBUTE_Species', right_on='organism_taxonomy_09species', how='left')
    df = pd.merge(df,df5,left_on= 'ATTRIBUTE_Genus', right_on='organism_taxonomy_08genus', how='left') 
    return df

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

def class_component(df1, df2):
    """ function to compute the class component based on the possible presence of new chemical classes 
    Args:
        df1 = reported_classes_df 
        df2 = sirius_classes_df
        Returns:
        None
    """
    if CC_component == False:
        print ('Similarit Class component not calculated')
    else:
        #merge the both tables
        df = pd.merge(df1,df2,on='filename', how='left').dropna()

        #get the difference between sets 

        df['New_in_species'] = df["CAN_classe"] - df["Chemical_class_reported_in_species"]  #check if the chemical classes from Sirius are reported in the species
        df['New_in_genus'] = df["New_in_species"] - df["Chemical_class_reported_in_genus"]  #check if the NEW chemical classes in the species are reported in the genus

        #Add the weight accordingly to the results 
        df['NS'] = df['New_in_species'].apply(is_empty)
        df['NG'] = df['New_in_genus'].apply(is_empty)
        
        #get the value of the CC 
        df['CC'] = df['NS'] + df['NG']
        
        df.drop('NS', axis=1, inplace=True)
        df.drop('NG', axis=1, inplace=True)
        return df