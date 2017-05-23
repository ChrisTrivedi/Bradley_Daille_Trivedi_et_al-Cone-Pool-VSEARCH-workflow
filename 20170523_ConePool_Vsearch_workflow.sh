#!/bin/bash

### Description ###
# This workflow uses a combination of QIIME and VSEARCH in order to create an OTU table that is suitable for input into multiple analysis pipelines (QIIME, R, etc.)

# This script assumes that you are in the directory that you are analyzing your data in. Paired end reads are located in Raw/ and subsequent files will be in seq/.

# While you can run this as a shell script, I find that running each command individually tends to work better.

# Also note the use of 'time' at the beginning of some of my scripts. This is purely for my own purposes to track how long individual scripts are taking to run.

### Dependent Software ###
# QIIME 1.9.1 at the time of this analysis
# VSEARCH (v2.3.4.osx_x86_64 used for this analysis), https://github.com/torognes/vsearch

### Workflow ###
# Starting with QIIME to join paired ends, extract barcodes, and split libraries.
# Join paired ends - overlap of 50
time join_paired_ends.py -f Raw/BS-S1_S1_L001_R1_001.fastq -r Raw/BS-S1_S1_L001_R2_001.fastq -o seq/Joined/ -j 50 -v

# Extract barcodes using the mapping file
time extract_barcodes.py -f seq/Joined/fastqjoin.join.fastq -m Parada_Amplicon_Run1_2016Map.txt -l 12 -o seq/Prepped/ -a -v

# Split libraries - separate out samples from fastq
time split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 0 -i Prepped/reads.fastq -b Prepped/barcodes.fastq -m Parada_Amplicon_Run1_2016Map.txt --barcode_type 12 -o SlOut/ -v

mkdir seq/VsearchOut

# Beginning of VSEARCH portion of the workflow
# get quality stats
vsearch -fastq_stats seq/SlOut/CPLMseqs.fastq --log seq/VsearchOut/seqs_stats.log

# remove low quality reads (trimming not required for paired-end data)
vsearch -fastq_filter seq/SlOut/CPLMseqs.fastq -fastaout seq/VsearchOut/seqs_filtered.fasta --fastq_maxee 0.5 --threads 4

# dereplicate seqs
vsearch -derep_fulllength seq/VsearchOut/seqs_filtered.fasta --output seq/VsearchOut/seqs_filtered_derep.fasta --sizeout --minuniquesize 2 --threads 4

# reference chimera check
vsearch -uchime_ref seq/VsearchOut/seqs_filtered_derep.fasta --db /macqiime/DB/gold.fasta --strand plus --nonchimeras seq/VsearchOut/seqs_filtered_derep_nochimeras.fasta --threads 4

# cluster OTUs @ 97%
vsearch -cluster_fast seq/VsearchOut/seqs_filtered_derep_nochimeras.fasta --centroids seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta --sizein --xsize --relabel OTU_ --id 0.97 --threads 4

# Make an otus folder
mkdir otus/

# Copy this file to use as RepSet at a later time
cp seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta otus/RepSet.fna

# map the original quality filtered reads back to OTUs
vsearch -usearch_global seq/VsearchOut/seqs_filtered.fasta --db seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta --strand plus --id 0.97 -uc  seq/VsearchOut/CPLM_otu_map.uc --threads 4

#Modify OTU table for input into QIIME
python /macqiime/bin/uc2otutab.py seq/VsearchOut/CPLM_otu_map.uc > seq/VsearchOut/seqs_filtered_derep_nochimeras_repset_OTU-table.txt

# convert to HDF5 biom type
biom convert --table-type="OTU table" -i seq/VsearchOut/seqs_filtered_derep_nochimeras_repset_OTU-table.txt --to-hdf5 -o otus/CPLM.biom

# Summarize BIOM table to check general stats
biom summarize-table -i otus/CPLM.biom -o CPLMSummary.txt

# Moving back into QIIME
# assign taxonomy
# You'll need to modify the "-t" flag to point to wherever your reference database is.  
# In this case I am using SILVA123 as my database of choice, but GreenGenes is the default for QIIME.
echo “Assigning Taxonomy”
time assign_taxonomy.py -t /macqiime/SILVA/taxonomy/taxonomy_all/97/taxonomy_7_levels.txt -r /macqiime/SILVA/rep_set/rep_set_all/97/97_otus.fasta -i otus/RepSet.fna -o otus/TaxonomyOut/

# add taxonomy to BIOM table
echo “Adding Metadata”
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/TaxonomyOut/RepSet_tax_assignments.txt -i otus/CPLM.biom -o otus/CPLM_otuTable.biom

#Computing Summaries of Taxa
echo "Computing Summaries"
biom summarize-table -i otus/CPLM_otuTable.biom -o CPLM_otuTable_Summary.txt

# Summarize the new BIOM table to check quality stats
biom summarize-table -i otus/CPLM_otuTable.biom --qualitative -o CPLM_otuTable_QualSummary.txt

# Align seqs (default: pynast - can also do in parallel)
time align_seqs.py -i otus/RepSet.fna -t /macqiime/SILVA/core_alignment/core_alignment_SILVA123.fasta -o otus/RepSet_Aligned/

# Filter alignment to certain quality
filter_alignment.py -i otus/RepSet_Aligned/RepSet_aligned.fasta -o otus/RepSet_Aligned/ -e 0.1

# Make tree file (.tre) for downstream use
make_phylogeny.py -i otus/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/RepSet_Aligned/rep_set.tre -l otus/RepSet_Aligned/tree_log.txt

### Done! ###
# The produced files can now be used for further analyses in QIIME or imported into an R package like 'phyloseq' (OTU table, mapping file, TRE file). Happy crunching!
