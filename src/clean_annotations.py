import pandas as pd
import numpy as np

## Clean annotations from GNPS
def cleaning_gnps_annotations(df1, df2, ppm_error, shared_peaks, cosine, ionisation_mode):
    """ function to clean and identify the valabe annotations from gnps

    Args:
        df1 = DB_results.tdv from gnps 
        df2 = Clusterinfosummary.tsv from gnps 
        ppm_error = min error in ppm to cosnider an annotation 'right'
        shared_peaks = min number of shared peaks between the ms2s spectra
        cosine = min cosine score between the ms2s spectra
        ionisation_mode = pos or neg according to experimental data
    
    Returns:
            pandas daraframe with the format of clusterinfosummary (all features from the network) + GNPS_INCHI_MF column (if not empty, the annotation for this feature was considered good) 

    """
    
    GNPS_annotations = df1
    print('==================')
    print('   Number of spectral library annotations in the job = '+str(df1.shape[0]))
    print('==================')

    #filter annotations not so good...
    GNPS_annotations = GNPS_annotations[GNPS_annotations.IonMode.str.lower().str.startswith(ionisation_mode, na=False)]
    GNPS_annotations = GNPS_annotations[GNPS_annotations.MZErrorPPM < ppm_error]
    GNPS_annotations = GNPS_annotations[GNPS_annotations.SpecCharge <= 1]
    GNPS_annotations = GNPS_annotations[GNPS_annotations.SharedPeaks > shared_peaks]
    GNPS_annotations = GNPS_annotations[GNPS_annotations.MQScore > cosine]
    
    # Entries that have salts or charged structure should be removed
    GNPS_annotations = GNPS_annotations[~GNPS_annotations.INCHI.str.contains('q:+|p+1|p+2|p-1', na=False)] 
    # Check that these are InChI
    GNPS_annotations = GNPS_annotations[GNPS_annotations.INCHI.str.upper().str.startswith(('INCHI','1S'), na=False)]
    # Filter the in source fragment spectral library match that are hard to match
    GNPS_annotations = GNPS_annotations[~GNPS_annotations.Adduct.str.contains("-C|i")]
    # Get the MF from the INCHI string
    GNPS_annotations['GNPS_INCHI_MF'] = GNPS_annotations.INCHI.str.split("/",expand=False).str[1]
    #Just to check there was no weird InChI
    GNPS_annotations = GNPS_annotations[GNPS_annotations.GNPS_INCHI_MF.str.upper().str.startswith('C', na=False)]   
    
    print('==================')
    print('   Number of spectral library annotations cleaned = '+str(GNPS_annotations.shape[0]))
    print('==================')

    #Merge the cleaned info with the information for the entire network 

    df =pd.merge(df2, GNPS_annotations[['#Scan#','GNPS_INCHI_MF']],left_on= 'cluster index', right_on='#Scan#', how='left')
    df.drop('#Scan#', axis=1, inplace=True)
    
    return df