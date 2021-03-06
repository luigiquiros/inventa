import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib
import memo_ms as memo 
import time

from sklearn.metrics import pairwise_distances
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest
from sklearn import preprocessing
from skbio.stats.ordination import pcoa
from skbio import OrdinationResults

#similarity component: 

def similarity_component(df, SC_component):
    """ function to compute the similarity component based on the MEMO matrix and machine learning unsupervised clustering methods 
    Args:
        df = memo matrix

    Returns:
        None
    """
    if SC_component == True:
        df1 = df.copy()
        df1.set_index('filename', inplace=True)
        df2 = df.copy()
        
        #specify the parameters of the individual classification algorithms
        clf = IsolationForest(n_estimators=100, 
                    max_samples='auto', 
                    contamination= 'auto', #0.15,
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
        df.to_csv('../data_out/SC_results.tsv', sep='\t')
        return df
    else:
        print('Similarity component not calculated')

def calculate_memo_matrix_ind_files(repository_path, spectra_suffix):
    
    # Generating memo matrix
    memo_unaligned = memo.MemoMatrix()

    start = time.process_time()
    memo_unaligned.memo_from_unaligned_samples(repository_path,  pattern_to_match =spectra_suffix, min_relative_intensity = 0.01,
                max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)
    print(f'Computing MEMO matrix from unaligned samples took: {time.process_time() - start} seconds')

    memo_unaligned.memo_matrix.index = memo_unaligned.memo_matrix.index.str.replace(spectra_suffix, "")

    memo_unaligned_filtered = memo_unaligned.filter(samples_pattern='01')
    memo_unaligned_filtered = memo_unaligned_filtered.filter(samples_pattern='12', max_occurence=0)
    
    df = memo_unaligned_filtered.memo_matrix
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'filename'}, inplace=True)
    return df 