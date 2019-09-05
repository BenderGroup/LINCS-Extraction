# LINCS Extraction
 Extraction of LINCS L1000 data

### Step 1: Create a folder called "data" in your working directory
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
--cell_line Specify cell line of interest for extraction
