<h1>Configurations and running</h1>

# Open the notebook 'inventa' 

Set the parameters according to your inputs and needs:

#### paths 

Each path corresponds to the files mentiones above. Just drop your files in the `/data` folder and change the names accordingly: 

```
    metadata_filename = '../data/Celastraceae_Set_metadata_pos.tsv'
    quantitative_data_filename = '../data/Celastraceae_pos_quant.csv'
    tima_results_filename = '../data/Celastraceae_pos_spectral_match_results_repond.tsv'
    vectorized_data_filename = '../data/Celastraceae_memomatrix.csv'
    canopus_npc_summary_filename = '../data/canopus_npc_summary.tsv'
    sirius_annotations_filename = '../data/canopus_npc_summary.tsv'

    #GNPS job id

    job_id= "yourjobidgoeshere"  #for example: job_id="4c919fcbc83d487493a487012afb920a"

```

#metadata

```
    species_column = 'yourspeciesnamecolumn'  #ATTRIBUTE_species
    genus_column =   'yourgenusnamecolumn'    #ATTRIBUTE_genus
    family_column =  'yourfamilynamecolumn'   #ATTRIBUTE_family
    organe_column =  'yourorganenamecolumn'   #ATTRIBUTE_organie
    filename_header = 'yourfilenamecolumn'    #filename
```

#### Parameters


##### For cleaning-up annotations from GNPS 

```
    ppm_error = 5                     # min error in ppm to consider an annotation valable
    shared_peaks = 10                 # min number of shared peaks between the MS2 experimental and MS2 fro, the database, to consider an annotation valable
    cosine = 0.7                      # min cosine score to consider an annotation valable
    ionisation_mode = 'pos'           # ionisation mode according to experimental conditions
```

##### Feature_component

```
    min_specificity = 90               # minimun feature specificity to consider
    only_feature_specificity = False   # True if annotations should be ignore and the FC should be calculated based on the features specificity. If False it will compute both The Sample specifity adn the FC

    #inputs to use: 

    isbd_annotations = True             # True: the tima_results_filename will be considered in the calculations
    sirius_annotations = True           #True: the sirius_annotations_filename will be considered in the calculations

    #cut-offs: 

    min_score_final = 0.0               #cut-off filter for considering an isdb annotation valable. You must be extremenly carefull with this parameter, '0.0' as default.
    min_ZodiacScore = 0.9               #cut-off filter for considering a sirius annotation valable. It is used in combination with min_ConfidenceScore.
    min_ConfidenceScore = 0.0           #cut-off filter for considering a sirius annotation valable. '0.0' as default.

    annotation_preference = 0          # Only Annotated nodes: '1' 
                                    # Only Not annotated: '0'
```

##### Literature_component

```
    LC_component = True                # LC will be calculated
    
    max_comp_reported = 40             # more than this value, the plant is considered no interesting LC =0
    min_comp_reported = 10             # less than this value, the plant is consireded very interesintg LC =1
                                    # a sample with x between both values gets a LC=0.5
```
##### Similarity_component

```
    SC_component = True                # SC will be calculated

```

##### Class_component

```
    
CC_component = True               # CC will be calculated
min_class_confidence = 0.8       #cut-off filter for considering a sirius class valable. It is used in combination with min_recurrence.
min_recurrence = 5               # minimum recurrence of a chemical class to consider it acceptable
```

##### specify the weight to modulate each component:

    # 1 means the value itself is taken into account. A 0.5 means onle half of the calculated value is taken into account

```
    w1 = 1   #for FC        
    w2 = 1   #for LC
    w3 = 1   #for SC
    w4 = 1   #for CC
```


## Once you set the parameter to your preference, Run the notebook!

# ENJOY! 


### [Back to Getting Started](getting-started.md)
### [Back to home page](index.md)

