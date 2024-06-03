Identify bacteriocin sequences present in clinical cohort of sequenced E. faecium isolates
This repository is split into 5 parts:
Part 1: Set up conda environment
Part 2: Snakefile to assemble and BLAST bacteriocin genes against isolates
Part 3: Filter BLAST hits 
Part 4: Cluster BLAST hits using CD-HIT
Part 5: Snakefile_variation to identify variants, low coverage regions, and insertion elements in hits to enterocin A and bacteriocin 43


# Part 1: Set up conda environment
Create an environment (here called snakemake2) and install snakemake and pyvcf
    
    mamba create -c conda-forge -c bioconda -n snakemake2 snakemake
    mamba install pyvcf3

Activate conda environment and load environment files

    conda activate snakemake2
    snakemake --use-conda --conda-create-envs-only -j96 

It may be necessary to add the conda bin to the path so conda can activate the environments

    export PATH="~/miniconda3/bin/:$PATH"


# Part 2: Snakefile
Edit config file to add shortread isolates under "short_read_samples" header, and samples for snippy core under "snippy_core_short_read_samples" header (see step 12 for more information)

    nano config/config.yaml


Make new directory "shortread" and add zipped shortread fastq files (should end in .fq.gz) 

    mkdir shortread

Do dry run of pipeline to make sure everything is working

    snakemake -n -r -s Snakefile -c24 --rerun-incomplete

Run pipeline

    sbatch ./scripts/snake_batch.sh 



#### Step 1: Trim shortreads 
Short read data is trimmed with trimmomatic . Parameters are indicated in the config.yaml. New data are deposted in a new directory 01_trimmed_shortread
#### Step 2: Assemble shortreads
Shor read data is assembled using Unicycler. Parameters are indicated in the config.yaml. New data are deposted in a new directory 02_assembly
#### Step 3: Index and annotate assemblies
Assemblies are annotated using Prokka. New data are deposited in a new directory 03_indexed_anno_reference.
#### Step 4: Identify MLST
Assemblies are annotated using mlst(https://github.com/tseemann/mlst). New data are deposited in a new directory 04_mlst.
#### Step 5: Make BLAST databases
Use BLAST command line to make a database from every isolate assembly. New data are deposited in a new directory 05_assembly_blast_databases.
#### Step 6: BLAST bacteriocins
BLAST a list of bacteriocin structural genes from Bagel4 against isolate assemblies. Bacteriocin genes are located in /blast_fastas/all_bacteriocin_classes.fa.  New data are deposited in a new directory 06_blast_all_bacteriocins.
#### Step 7: BLAST enterocin A
BLAST the genes associated with enterocin A production against isolate assemblies. Genes are located in /blast_fastas/enterocinA.fa.  New data are deposited in a new directory 07_blast_enterocinA.
#### Step 8: BLAST bac43
BLAST the genes associated with bacteriocin 43 production against isolate assemblies. Genes are located in /blast_fastas/bac43.fa.  New data are deposited in a new directory 08_blast_bac43.
#### Step 9: BLAST ent53B
BLAST the genes associated with enterocin-5-3B production against isolate assemblies. Genes are located in /blast_fastas/ent53B.fa.  New data are deposited in a new directory 09_blast_ent53B.
#### Step 10: BLAST van genes
BLAST the genes associated with vanomycin resistance (vanA, vanB, vanM, vanD). Genes are located in /blast_fastas/vanABMD.fa.  New data are deposited in a new directory 10_blast_vanABMD.
#### Step 11: Identify variants
Identify variants using snippy  Parameters are indicated in the config.yaml. Reference is located in snippy_reference/BL00198-1.gbk. New data are deposited in a new directory 11_snippy.
#### Step 12: Create a core genome
Use snippy_core to create a core genome. Parameters are indicated in the config.yaml. Reference is located in snippy_reference/BL00198-1.gbk.
New data are deposited in a new directory 12_core_genome.
If snippy_core rule fails with 
    Warning: No SNPs were detected so there is nothing to output.
Use text file produced by failed snippy_core rule (should be ouputed to 12_core_genome/{reference}.txt to remove samples with less than 250,000bp aligned. Create a new list of isoaltes with greater than 250,000bp aligned and add to config/config.yaml under the "snippy_core_short_read_samples" header.
#### Step 13: Generate a phylogenetic tree
Use gubbins to create a phylogenic tree. New data are deposited in a new directory 13_gubbins.

# Part 3: Filter BLAST output 
Use R to run 

    analysis/filter_blast_output.Rmd

This will filter BLAST outputs by e-value and remove multiple query hits to the same contig region.

# Part 4: Cluster BLAST hits
Set up conda environment for cdhit and install
    
    conda install -c bioconda cd-hit

Use R to run first part of 

    analysis/process_cluster_blast_output.Rmd

This will merge blast output with list of isolates, pids and mlsts, filter hits to those with greater than 60% identity, and create a fasta file with hits for clustering 
(first_16to21_besthit_sequences.fa)

Activate conda environment for cdhit and run

    conda activate cdhit
    cd-hit -i first_16to21_besthit_sequences.fa -o cdhit_clustered_hits90.fa -c 0.9 -G 1  -l 10 -g 1 -aL 0.75 -sc 1

Go back to R and run second part of 

    analysis/process_cluster_blast_output.Rmd

This will assign clusters to hits, reassign clusters of split structural genes for bacteriocins, reassign clusters numbers (for clusters that were removed so that they are sequential), rename queries, and output clustered hits for analysis - clustered_hits.csv











