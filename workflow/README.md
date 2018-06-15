VARSCOT workflow
---

**Scripts**

- processDataForModel.R: prcess GUIDE-seq data and mapping results from bidirectional index for inactive off-targets in order to create 10 datasets for training. Data is filtered and in the end inactive off-targets are sampled and final datasets are constructed.
- classificationModel.R: builds feature matrix for all 10 training datasets, performs feature selection, builds models and compares classifiers based on cross validation. Best performing model is selected as the final model.
- siteseqBiochemicalValidation.R: performs validation of the model on the SITE-seq data
- siteseqPipelineVomparison.R: perfoms pipelien comparisons using the SITE-seq data
- evalFunctions.R: various functions used by multiple notebooks.

**Folders**

- data-objects: Contains RData objects for and from the notebooks
- guideseq-data: Contains data from the GUIDE-seq dataset
- siteseq-data: Contains data from the SITE-seq dataset
- pipeline-comparison: Contains data for the pipeline comparison

