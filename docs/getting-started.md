<h1>Getting started</h1>


To run the repository the minimun inputs needed are:

- GNPS job
- Metadata table (the same used for the GNPS job)
- Quantification table from MZmine or another similar software

This will run at least the Feature component (FC).

Optionally, the following inputs can be added: 

- Annotations files (ISDB and/or SIRIUS)
- Canopus file (chemical ontology)
- vectorized dissimilarity matrix file (MEMO)


# Where to start? 


# Check the inputs 


#### 1.1 Metadata table:
The standard format from GNPS is prefered:

`metadata_filename`: it uses the [GNPS format](https://docs.google.com/spreadsheets/d/1pSrqOdmMVBhVGpxIZeglToxihymTuaR4_sqTbLBlgOA/edit#gid=0).

While creating the 'metadata' there are some MANDATORY columns:

Taxonomy: the species, genus and family are neeeded if the LC and CC want to be computated. The taxonomy should be cleaned to uptoday recognized names, you can use the [Open Tree of Life](https://opentree.readthedocs.io/en/latest/readme.html).

The headers for each one could follow the GNPS format or the user's preferences, how ever the following parameter need to be indicated:

```
        species_column = 'yourspeciesnamecolumn'  #ATTRIBUTE_species
        genus_column =   'yourgenusnamecolumn'    #ATTRIBUTE_genus
        family_column =  'yourfamilynamecolumn'   #ATTRIBUTE_family
        organe_column =  'yourorganenamecolumn'   #ATTRIBUTE_organie
        filename_header = 'yourfilenamecolumn'    #filename
```

The `organe _colum` should be specified if you have diferent parts (or solvents) from the same species. If you prefer to use only the filename as identifier for the resuts, it can be specified directly in the notebook.

#### 1.2 Feature quantitative table:

`quantitative_data_filename`: MZmine output format using only the 'Peak area', 'row m/z' and 'row retention time' columns.  

- if you prefer 'Peak Height', go to `src/process_data.py` and change it inside the function quant_table(). ONLY ONE of the columns is considered at the time, 'Peak height' or 'Peak area', if you want to consider both they must be done one at a time.  

- if you did export any other column, like identities, etc,  please remove manually or add the corresponding lines in the funcion quant_table():

```
        df.drop('name of the colum', axis=1, inplace=True)
```

- Usualy there are columns with the header 'Unkown: number' at the very end of the quantitative table, the script takes care of theses columns, you do not need to erase them

#### 1.3 In silico annotation usign timaR:

`tima_results_filename`: [timaR](https://taxonomicallyinformedannotation.github.io/tima-r/index.html) reponderated output format.

- for performing in silico annotations and taxonomically informed reponderation.

#### 1.4 Chemical taxonomy results:

`canopus_npc_summary_filename`: Sirius CANOPUS recomputated output format.

This output needs an additional step after runnign sirius, please follow the next instructions:

- install  [Sirius](https://bio.informatik.uni-jena.de/software/sirius/) and run it in your set.
- clone [Canopus treemap](https://github.com/kaibioinfo/canopus_treemap)
- Recompute your project space from Sirius using the following code:

``` 
        from canopus import Canopus
        C = Canopus(sirius="sirius_projectspace")
        C.npcSummary().to_csv("npc_summary.tsv")
```

- the output `canopus_npc_summary.tsv` corresponds to the file nedded for running Inventa

- given that the [Lotus Dabase](https://lotus.naturalproducts.net/) uses the NPClassifyre ontology and Sirius uses the Classifyre ontology, performing this step is absolutley necesary for a proper comparison of the propsed chemical classes.

#### 1.5 Annotations with Sirius: 

`sirius_annotations_filename`: Sirius annotations output format. Containing Zodiac and Cosmic.

- this file should correspond to `compound_identification.tsv`
- make sure the Zodiac and Cosmic scores are present ().

#### 1.6 Memo dissimilarity matrix:

`vectorized_data_filename`: [MEMO](https://github.com/mandelbrot-project/memo) package output format.

Examples of all these input could be found in [`/format_examples`](https://github.com/luigiquiros/inventa/tree/main/format_examples)



## [Continue to Configurations and running](configuration-options.md) 

### [Back to home page](index.md)