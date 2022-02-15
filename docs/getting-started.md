<h1>Getting started</h1>


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


# Check the inputs 


#### 1.1 Metadata table:
The standard format from GNPS is prefered:

`metadata_filename`: GNPS format ([https://docs.google.com/spreadsheets/d/1pSrqOdmMVBhVGpxIZeglToxihymTuaR4_sqTbLBlgOA/edit#gid=0](https://docs.google.com/spreadsheets/d/1pSrqOdmMVBhVGpxIZeglToxihymTuaR4_sqTbLBlgOA/edit#gid=0)).

While creating the 'metadata' there are some MANDATORY headers:

- `ATTRIBUTE_Species` : The species should be cleaned to uptoday recognized names, you can use the Open Tree of Life to do so (https://opentree.readthedocs.io/en/latest/readme.html).
- `ATTRIBUTE_Organe`  : This column correpond to the part of the plant or organism.

#### 1.2 Feature quantitative table:

`quantitative_data_filename`: MZmine output format using only the 'Peak area', 'row m/z' and 'row retention time' columns.  

- if you prefer 'Peak Height', go to `src/inventa.py`and change it inside the function quand_table(). ONLY ONE of the columns is considered at the time, 'Peak height' or 'Peak area', if you want to consider both they must be done one at a time.  

- if you did export any other column, like identities, etc,  please remove manually or add the corresponding lines in the funcion quand_table():
```
 df.drop('name of the colum', axis=1, inplace=True)
 ```

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

[Back to home page](index.md)