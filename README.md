# LINCS Extraction
 Extraction of LINCS L1000 data.
 
 [A Next Generation Connectivity Map: L1000 platform and the first 1,000,000 profiles](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5990023/)
 
 Repo image from [The Library of Integrated Network-Based Cellular Signatures NIH Program: System-Level Cataloging of Human Cells Response to Perturbations](https://www.ncbi.nlm.nih.gov/pubmed/29199020)
 
 Authors: Maria-Anna Trapotsi (mat64@cam.ac.uk) and Layla Hosseini-Gerami (lh605@cam.ac.uk)

Before you start, make sure you have the following packages installed. We recommend you create a Conda environment for this purpose.

* pandas
* numpy
* cmapPy
* h5py
* scipy
* joblib

### Step 1: Create folders
Create folders called "data", "Consensus_Signatures", "TAS_Scores", "Compounds_With_No_Replicates", "All_Replicates" in your working directory

 ```
Working Directory/
├── data/
│   ├── GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
│   ├── GSE92742_Broad_LINCS_gene_info.txt.gz
│   ├── GSE92742_Broad_LINCS_sig_info.txt.gz 
│   ├── GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
│   └── GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
├── Consensus_Signatures/
├── TAS_Scores/
├── Compounds_With_No_Replicates/
├── All_Replicates/
├── automation.sh
└── Extract_And_TAS.py
```

### Step 2: Extract LINCS L1000 files from Phase 1 and Phase 2 from Gene Expression Omnibus (GEO) and save them in "data" folder
For Phase 1 visit:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
and download the following files:
GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
GSE92742_Broad_LINCS_sig_info.txt.gz
GSE92742_Broad_LINCS_gene_info.txt.gz

For Phase 2 visit:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138
and download the following files:
GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz

### Step 3: Run the Extract_And_TAS.py

### Extract_And_TAS.py options:

```
Options:
  -h, --help            show this help message and exit
  --cell_line=CELL_LINE
                        Cell line for extraction (default MCF7)
  --pert_time=PERTURBATION_TIME
                        Perturbation time for extraction (hours) (default 24)
  --pert_dose=PERTURBATION_DOSE
                        Perturbation dose for extraction (uM) (default 10)
  --ncores=CORES        Number of cores (default 1)
 ```
 
I.e. if you want to extract all data from VCAP measured in 24 hours and 5 uM, using 5 cores of compute power you would type:

```
$ python Extract_And_Tas.py --cell_line=VCAP --pert_time=24 pert_dose=5 --ncores=5
```

Please note you may have to change the filepaths in the code to match with the filenames you have extracted in data/ folder, as filenames are periodically updated.
 
**This script will extract LANDMARK GENES ONLY but you can modify a couple of lines in the code to extract best inferred (BING) or all genes also. Please get in touch if you need help with this.**
  
automation.sh will enable you to extract data in multiple conditions, by providing a few parameters which will be looped over in every combination - you may edit this file to keep only cell lines etc. you care about.

```
$ chmod +x automation.sh # make it executable
$ ./automation.sh
```

 ### Extracted Data
 1. Compounds with replicates: 
Consensus signatures are saved in the full matrix of all compounds in conditions specified (e.g A375_24 h_10uM.txt) where headers are genes and rows are compounds, in the "Consensus_Signatures" folder.

Consensus signatures are computed for biological replicates in the same way that the original paper details for technical replicates.

Individual replicates are also extracted; these are saved as a matrix per compound (e.g. compoundname_A375_24 h_10uM.txt) where rows are genes and headers are compound replicate ID, in the "All_Replicates" folder.

2. Compounds with no replicates:
Individual signatures are saved in the full matrix of all compounds in conditions specified (as above) but in the "Compounds_With_No_Replicates" folder.

3. TAS scores for all compounds with replicates
TAS score + SS and CC scores for all compounds with replicates in conditions specified (e.g. A375_24 h_10uM_TAS.txt) in TAS_Scores folder.

TAS scores are computed for biological replicates in the same way that the original paper details for technical replicates.

**!!!If you just want a matrix of gene expression data then you will need to combine the matrix in Consensus_Signatures and Compounds_With_No_Replicates for your condition of interest (simple Python concat function)!!!**

### Future functionality
1. Ability to provide more than one cell line, time point, dose etc. as an option (for now automation.sh script allows iteration, but each extraction is a separate task i.e. not as a loop within the Python script)
2. Command line option to choose gene subset (LM, BING, All - for now just LM but this can manually be edited in the Python script)
3. Command line option for TAS score calculation or not



                           
                             
                                           
                           
