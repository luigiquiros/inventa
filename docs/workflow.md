# Workflow 

To run the repository the minimun inputs needed are:

- GNPS job
- Metadata table (the same used for the FBMN)
- Quantification table

This will run at least the Feature component (FC).

Optionally, the following inputs can be added: 

- Annotations files (ISDB and/or SIRIUS)
- Canopus file (chemical ontology)
- vectorized dissimilarity matrix file (MEMO)


# Where to start? 


## 1. Check the inputs 


#### 1.1 Metadata table:
The standard format from GNPS is prefered:

`metadata_filename`: GNPS format ([https://docs.google.com/spreadsheets/d/1pSrqOdmMVBhVGpxIZeglToxihymTuaR4_sqTbLBlgOA/edit#gid=0](https://docs.google.com/spreadsheets/d/1pSrqOdmMVBhVGpxIZeglToxihymTuaR4_sqTbLBlgOA/edit#gid=0)).

While creating the 'metadata' there are some MANDATORY headers:

- `ATTRIBUTE_Species` : The species should be cleaned to uptoday recognized names, you can use the Open Tree of Life to do so (https://opentree.readthedocs.io/en/latest/readme.html).
- `ATTRIBUTE_Organe`  : This column correpond to the part of the plant or organism.

#### 1.2 Feature quantitative table:

`quantitative_data_filename`: MZmine output format using only the 'Peak area', 'row m/z' and 'row retention time' columns.  

- if you prefer 'Peak Height', go to `src/inventa.py`and change it inside the function quand_table(). ONLY ONE of the columns is considered at the time, 'Peak height' or 'Peak area', if you want to consider both they must be done one at a time.  

- if you did export any other column, like identities, etc,  please remove manually or add the corresponding lines in the funcion quand_table(), `df.drop('name of the colum', axis=1, inplace=True)`.
- Usualy there are columns with the header 'Unkown: number' at the very end of the quantitative table, the scrip takes care of theses columns, you do not need to erase them

#### 1.3 In silico annotation usign timaR:

`tima_results_filename`: timaR reponderated output format.

- for performing in silico annotations and taxonomically informed reponderation, please follow https://taxonomicallyinformedannotation.github.io/tima-r/index.html

#### 1.4 Chemical ontology results:

`canopus_npc_summary_filename`: Sirius CANOPUS recomputated output format.

This output needs an additional step after runnign sirius, please follow the next instructions:

- if you don't have Sirius, please install it from here (https://bio.informatik.uni-jena.de/software/sirius/), and run it in your set. 
- clone the following repository https://github.com/kaibioinfo/canopus_treemap
- Recompute your project space from Sirius using the following code:

``` 
        from canopus import Canopus
        C = Canopus(sirius="sirius_projectspace")
        C.npcSummary().to_csv("npc_summary.tsv")
```

- the output `canopus_npc_summary.tsv` corresponds to the file nedded for running Inventa

- given that the Lotus Dabase (https://lotus.naturalproducts.net/) uses the NPClassifyre ontology and Sirius uses the Classifyre ontology, performing this step is absolutley necesary for a proper comparison of the propsed chemical classes.

#### 1.5 Annotations with Sirius: 

`sirius_annotations_filename`: Sirius annotations output format. Containing Zodiac and Cosmic results (https://bio.informatik.uni-jena.de/software/sirius/).

- this file should correspond to `compound_identification.tsv`
- make sure the Zodiac and Cosmic scores are present ().

#### 1.6 Memo dissimilarity matrix:

`vectorized_data_filename`: MEMO package format (https://github.com/mandelbrot-project/memo).

[Examples of all these input could be found in `/format_examples`]



## 2. Open the notebook 'inventa' 

Set the parameters according to your inputs and needs:

#### 2.1 paths 

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

#### 2.2 Parameters


#### 2.2.1 For cleaning-up annotations from GNPS 

```
ppm_error = 5                     # min error in ppm to consider an annotation valable
shared_peaks = 10                 # min number of shared peaks between the MS2 experimental and MS2 fro, the database, to consider an annotation valable
cosine = 0.7                      # min cosine score to consider an annotation valable
ionisation_mode = 'pos'           # ionisation mode according to experimental conditions
```

#### 2.2.1 Feature_component

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

#### 2.2.2 Literature_component

```
LC_component = True                # LC will be calculated
max_comp_reported = 40             # more than this value, the plant is considered no interesting LC =0
min_comp_reported = 10             # less than this value, the plant is consireded very interesintg LC =1
                                   # a sample with x between both values gets a LC=0.5
```
#### 2.2.3 Similarity_component

```
SC_component = True                # SC will be calculated

```

#### 2.2.4 Class_component

```
CC_component = True                # CC will be calculated
min_recurrence = 10                # minimun recurrence of a chemical class to considered it valable.
```

#### 2.2.5 specify the weight to modulate each component 

```
w1 = 1           # 1 means the value itself is taken into account. A 0.5 means onle half of the calculated value is taken into account
w2 = 1
w3 = 1
w4 = 1
```


## 3. Run the notebook 'inventa' and ENJOY! 