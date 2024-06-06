Identify bacteriocin sequences present in clinical cohort of sequenced E. faecium isolates  
This repository is split into 6 parts:  
**Part 1**: Set up conda environment  
**Part 2**: Snakefile to assemble and BLAST bacteriocin genes against isolates  
**Part 3**: Filter BLAST hits   
**Part 4**: Cluster BLAST hits using CD-HIT  
**Part 5**: Snakefile_variation to identify variants, low coverage regions, and insertion   elements in hits to enterocin A and bacteriocin 43  
**Part 6**: Filter and merge variations outputed by part 5


# Part 1: Set up conda environment
Create an environment (here called snakemake2) and install snakemake and pyvcf
    
    mamba create -c conda-forge -c bioconda -n snakemake2 snakemake
    mamba install pyvcf3

Install plugin to use cluster

    mamba install snakemake-executor-plugin-cluster-generic


# Part 2: Snakefile
Edit config file to add shortread isolates under "short_read_samples" header, and samples for snippy core under "snippy_core_short_read_samples" header (see step 12 for more information)

    nano config/config.yaml


Make new directory "shortread" and add zipped shortread fastq files (should end in .fq.gz) 

    mkdir shortread

If doing a hybrid assembly, also make a directory "longread" and add zipped longread fastq files. Add isolates to the config file under the header "long_read_samples".

Activate conda environment and load environment files

    conda activate snakemake2
    snakemake --use-conda --conda-create-envs-only -j96 

It may be necessary to add the conda bin to the path so conda can activate the environments

    export PATH="~/miniconda3/bin/:$PATH"

Do dry run of pipeline to make sure everything is working

    snakemake -n -r -s Snakefile -c24 --rerun-incomplete

Run pipeline

    sbatch ./scripts/snake_batch.sh 



#### Step 1: Trim shortreads 
Short read data is trimmed with trimmomatic . Parameters are indicated in the config.yaml. New data are deposted in a new directory 01_trimmed_shortread
#### Step 1a: Trim longreads 
(For hybrid assembly) Long read data is trimmed with porechop . Parameters are indicated in the config.yaml. New data are deposted in a new directory 01_trimmed_longread
#### Step 2: Assemble shortreads
Short read data is assembled using Unicycler. Parameters are indicated in the config.yaml. New data are deposted in a new directory 02_assembly
#### Step 2a: Hyrbid assembly
(For hybrid assembly) Short read and long read data is assembled using Unicycler. Parameters are indicated in the config.yaml. New data are deposted in a new directory 02_hybrid_assembly
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

Use text file produced by failed snippy_core rule (should be ouputed to 12_core_genome/{reference}.txt) to remove samples with less than 250,000bp aligned. Create a new list of isoaltes with greater than 250,000bp aligned and add to config/config.yaml under the "snippy_core_short_read_samples" header.

    scripts/determine_snippy_core_isolates.Rmd
    
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


# Part 5: Snakefile_variation
Activate conda environment and load environment files

    conda activate snakemake2
    snakemake -s Snakefile_variation --use-conda --conda-create-envs-only -j96

It may be necessary to add the conda bin to the path so conda can activate the environments

    export PATH="~/miniconda3/bin/:$PATH"
Do dry run of pipeline to make sure everything is working

    snakemake -n -s Snakefile_variation -c24 

Run pipeline

    sbatch ./scripts/snake_batch_variation.sh 

### Define variation in enterocin A and bactericoin 43
Each of these steps (except where indicated) were repeated for both enterocin A and bacteriocin 43. Different internal references were used for each bacteriocin.   
Bacteriocin 43: bwa_references/BL02040-1_assembly.fasta   
Enterocin A: bwa_references/PR46485-1-C1_assembly.fasta   

#### Step 1: Index reference
Use BWA to index reference files
#### Step 2: Map shortreads
Use BWA to map shortread files to regions of interest for bacteriocins.   
Bacteriocin 43: Contig 7
Enterocin A: Contig 2,  246177 - 256268   
#### Step 3: Identify and annotate variants 
Use freebayes to identify variants and snpEFF to annotate variants.
#### Step 4: Identify and annotate insertion elements
Use panISa to identify potential insertion elements. Then run ISFinder through script adopted from panISa to annotate variants (scripts/ISFinder_search.py) 
#### Step 5: Identify areas of low coverage (deletions)
Use samtools to identify areas of low coverage as potential gene deletions
### Step 6: Create IGV report
Use IGV to output a read pileup map over bacteriocin regions of interest 
### Step 7: Get coverage (bacterioin 43 only)
Use samtools to calculate coverage across the entire bac43 plasmid and across the bac43 genes 


# Part 6: Filter and merge variations
After manually inspecting IGV read pileups and confirming variants, use R to create a master list of variants, insertion elements, and deletions. 

    scripts/merge_variations_enterocinA.Rmd
    scripts/merge_variations_bac43.Rmd

For bacteriocin 43, isolates and variants were evaluated using three rules:   
1. Filter out isolates that have sparse coverage to plasmid (45% or greater mean depth across the plasmid and 55% or greater coverage)   
2. Split out positions that have split hits to the bac43 plasmid and other regions of the genome (DPtoRO < 0.2  (i.e. reads calling the reference allele make up less than 20% of the total reads mapped) and RO < 100 (i.e. less than 100 reads mapping to reference allele))   
3. Identify true variants to the bacteriocin 43 plasmid (percent identity greater than 99.5% unless bacteriocin 43 was identified in the BLAST search)   

    












