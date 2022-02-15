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

#other: 

only_ms2_annotations = False       # False to considere both, MS1 & MS2 annotations, False will only considerer MS2 annotations
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
CC_component = True                # CC will be calculated
min_recurrence = 10                # minimun recurrence of a chemical class to considered it valable.
```

##### specify the weight to modulate each component:

    # 1 means the value itself is taken into account. A 0.5 means onle half of the calculated value is taken into account
```
w1 = 1   for FC        
w2 = 1   for LC
w3 = 1   for SC
w4 = 1   for CC
```


# Run the notebook 'inventa' and ENJOY! 

[Back to home page](index.md)
[Back to Getting Started](getting-started.md)
