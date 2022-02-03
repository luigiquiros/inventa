# INVENTA


## Description 

This is a workflow to explore the potential of a set of samples to contain not previously reported compounds. It is based on the results from data treatment usign MZmine, spectral organization through Molecular Netwroking, in-silico dereplication and prediction. 
It is composed of 4 independend components: 

The **feature component (FC)** is a ratio that considers compounds with a specificity higher than, for instance, 90% per sample (specified by the user), and without putative annotation. Results: 1 Column with  FC ratio and 1 column with the sample specificity (ratio of peaks higher thant the specificied % wth or without annotations).

The **literature component (LC)** is a score based on the number of compounds reported in the literature for the taxon. The output includes one column with the LC score and at least two additional columns of metadata containing the number of reported compounds at the species and genus levels. 

The **class component (CC)** is a score based on the presence of possible new chemical classes in the taxon, not previously reported before. The CC will be considered an integer ‘1’ if there are new chemical classes at the species level, and an additional ‘1’ if those new chemical classes are not present in the genus either.

The similarity component (SC) is a score based on the spectral similarity of the sample within the set. Multiple outlier detection machine learning algorithms are implemented to spot the dissimilar samples. A weight of ‘1’ is given to the sample considered anomalies in at least one detection method.

### Running INVENTA with Binder

Binder allows to run INVENTA on the cloud with a Binder instance, which is really convenient but you need to save the parameters and results locally as these instance are shutting down after 15 min of inactivity.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/luigiquiros/INVENTA/main?labpath=notebooks%2INVENTA_v7.ipynb)


### Running INVENTA locally:

#### Install the conda environment

First make sure to have [anaconda](https://www.anaconda.com/products/individual) installed.

#### Clone the repo locally
First clone the repository using git clone in command line. You may need to install the git package (see [here](https://www.atlassian.com/git/tutorials/install-git):
```
git clone https://github.com/luigiquiros/INVENTA.git
```

Create a new conda environment to avoid clashes:
```
conda env create -n inventa -f inventa/environment.yml
```
Use this one to activate the environment: 
```
conda activate inventa
```

If you need to update the environment run, make sure to run the following to keep the depencies versions
```
conda env update --file environment.yml
```

If you have an error, try installing `scikit-bio` from `conda-forge` before installing the package with pip:
```
conda install -c conda-forge scikit-bio
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ memo-ms
```



### The format for imput tables is critical! please read carefully the following lines to asure you have the right format
#### Metadata table format
Use the format of the metadata tables specified:
    `metadata_filename`: GNPS format ([https://docs.google.com/spreadsheets/d/1pSrqOdmMVBhVGpxIZeglToxihymTuaR4_sqTbLBlgOA/edit#gid=0](https://docs.google.com/spreadsheets/d/1pSrqOdmMVBhVGpxIZeglToxihymTuaR4_sqTbLBlgOA/edit#gid=0)).
    - the `ATTRIBUTE_Species` and `ATTRIBUTE_Organe` headers are mandatory. The species should be cleaned to uptoday recognized names.
    - the `ATTRIBUTE_Sppart` is generated in the notebook from the ATTRIBUTE_Species and ATTRIBUTE_organe colums, if you already have this column in your metadata be sure the header match properly and ignore the line in the data preparation section. 

#### Feature quantitative table
  - `quantitative_data_filename` = MZmine output format using only the 'Peak area', 'row m/z' and 'row retention time' columns.  
-if you prefer 'Peak Height', go to INVENTA > src > inventa.py and change it inside the function quand_table(). ONLY ONE of the columns is considered at the time, 'Peak height' or 'Peak area', if you want to consider both they must be done one at a time.  
-if you did export any other column, like identities, etc,  please remove manually or add the corresponding lines in the funcion quand_table(), `df.drop('name of the colum', axis=1, inplace=True)`  

#### Annotation tables 
Normally there are columns with the header 'Unkown: number' at the very end of the quantitative table, the scrip takes care of theses columns, you do not need to erase them.  
    `clusterinfosummary` : GNPS format as downloaded from the job.
    `isdb_results_filename` : format from the metadata_annotator, JLW lab.
    `vectorized_data_filename` : MEMO package format.
    `sirius_results_filename` : CANOPUS/SIRIUS format. npc_summary_network

There are some parameters that need to be fixed by the user before launching the job. 
GO TO INVENTA > src > inventa.py and cange accordingly: 

### INVENTA parameters

#### Feature component

        FC_component = True                          #FC will be calculated
        min_specificity = 90                         #minimun feature specificity to consider
        only_feature_specificity = False             #True if annotations should be ignore and the FC should be calculated based on the features specificity. If False it will compute both The Sample specifity adn the FC
        only_gnps_annotations = False                #only the annotations from gpns will be considered 
        only_ms2_annotations = False                 #False to considere both, MS1 & MS2 annotations, False will only considerer MS2 annotations
        annotation_preference = 0                     #Only Annotated nodes: '1' /  Only Not annotated: '0'

#### Literature component 

        LC_component = True                         #LC will be calculated
        max_comp_reported = 40                      #more than this value, the plant is considered no interesting LC =0
        min_comp_reported = 10                      #less than this value, the plant is consireded very interesintg LC =1, a sample with x between both values gets a LC=0.5
        family_compounds = False                    #True is the nomber of reported in the family should be retreived

### Class component+
`CC_component = True  #CC will be calculated`

### Similarity component
`SC_component = True  #SC will be calculated`

Drop your files in the data folder and change the names in the notebook to march them: 

### Input filenames: drag them in the data folder
```
    metadata_filename = '../data/Celastraceae_Set_metadata_pos.tsv'
    quantitative_data_filename = '../data/Celastraceae_pos_quant.csv'
    clusterinfosummary = '../data/198a574e172443ad84f43e739ceea2c8.tsv'
    isdb_results_filename = '../data/Celastraceae_pos_spectral_match_results_repond.tsv'
    vectorized_data_filename = '../data/Celastraceae_memomatrix.csv'
    sirius_results_filename = '../data/canopus_npc_summary_CANOPUS_network.txt'
```
### ENJOY!!! 
