# create_concensus_genome
Script to create concensus genome

## Step 1: get_sequences.py
Query accession from UCSC API

## Step 2: check_anchor.py
Use data from ncbi to get which sequence is anchor or not
 
## Step 3: get_coord_onhg38.py
Search position which is anchor on hg38 to extract major 

## Step 4: extract_major.py 
Extract major variants from major file and change name and pos to corresponding alt contig

## Step 5: 
Use bcftools consensus to replace the major variants in GRCh38 genome