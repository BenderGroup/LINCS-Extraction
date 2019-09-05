# LINCS Extraction
 Extraction of LINCS L1000 data

### Step 1: Create folders
Create folders called "data", "Consensus_Signatures", "TAS_Scores", "Compounds_With_No_Replicates", "All_Replicates" in your working directory

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
                        Perturbation time for extraction (hours) (default 24h)
  --pert_dose=PERTURBATION_DOSE
                        Perturbation dose for extraction (uM) (default 10)
  --ncores=CORES        Number of cores (default 1)
 ```
  
 ### Extracted Data
 1. Compounds with replicates: 
Consensus signatures are saved in the full matrix of all compounds in conditions specified (e.g A375_24 h_10uM.txt) where headers are genes and rows are compounds, in the "Consensus_Signatures" folder.
Individual replicates are also extracted; these are saved as a matrix per compound (e.g. compoundname_A375_24 h_10uM.txt) where rows are genes and headers are compound replicate ID, in the "All_Replicates" folder.

2. Compounds with no replicates:
Individual signatures are saved in the full matrix of all compounds in conditions specified (as above) but in the "Compounds_With_No_Replicates" folder.

3. TAS scores for all compounds with replicates
TAS score + SS and CC scores for all compounds with replicates in conditions specified (e.g. A375_24 h_10uM_TAS.txt) in TAS_Scores folder.

### Future functionality
1. Ability to provide more than one cell line, time point, dose etc. as an option (for now automation.sh script allows iteration, but each extraction is a separate task i.e. not as a loop within the Python script)
2. Command line option to choose gene subset (LM, BING, All - for now just LM but this can manually be edited in the Python script)
3. Command line option for TAS score calculation



                           
                                    |LINCS|
                                    /  |  \
                                   /   |   \
                                  /    |    \
                                 /     |     \
    |Compounds_With_No_Replicates|   |DATA|  |Consensus_Signatures|
                                    
                           
