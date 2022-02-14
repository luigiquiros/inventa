import pandas as pd
import numpy as np
import zipfile
import os
import scipy as sp
import matplotlib.pyplot as plt
import plotly.express as px
import zipfile
import pathlib

def get_gnps_annotations(df):
    #retrive the clusterinfosummary file from the gnps jod downloaded before
    
    for filepath in pathlib.Path("../data/all_annotations/clusterinfo_summary/").glob('**/*'):
        filepath.absolute()
    
    df_network = pd.read_csv(filepath.absolute(),
                                sep='\t', 
                                usecols =['cluster index', 'componentindex'])
    df = pd.merge(df_network[['cluster index', 'componentindex']], df,left_on= 'cluster index', right_on='#Scan#', how='left')
    
    return df

def get_isdb_annotations(path_isdb, isdb_annotations): 
    
    if isdb_annotations == True:
        df = pd.read_csv(path_isdb,
                                sep='\t', 
                                usecols =['feature_id','molecular_formula','score_final','score_initialNormalized'], 
                                low_memory=False)
        return df 
    else: 
        print('The isdb annotations output will be not used')

def get_sirius_annotations(path_sirius, sirius_annotations): 

    if sirius_annotations == True:
        df = pd.read_csv(path_sirius,
                                sep='\t', 
                                usecols =['id','ConfidenceScore','ZodiacScore'], 
                                low_memory=False)
        return df 
    else: 
        print('The sirius annotations will be not used')

def get_canopus_pred_classes(path_canopus, CC_component): 

    if CC_component == True:
        df = pd.read_csv(path_canopus, sep='\t')
        return df 
    else: 
        print('The canopus classes will be not used')

