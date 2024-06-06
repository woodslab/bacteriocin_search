###Snakefile for sample assembly and bacteriocin blast
configfile:"config/config.yaml"

SAMP=config["short_read_samples"]
SAMP_CORE=config["snippy_core_short_read_samples"]
LONG_SAMP=config["long_read_samples"]
INDEX_FILE = ["amb", "ann", "bwt", "fai", "pac", "sa"]
REFERENCE=config["reference_genome"] #for snippy

rule all:
    input: 
        expand("01_trimmed_shortread/{sample}_01.fq.gz", sample=SAMP),
        #expand("01_trimmed_longread/{sample}_longread.fq", sample=LONG_SAMP),
        expand("02_assembly/{assembly_sample}/{assembly_sample}_assembly.fasta", assembly_sample=SAMP),
        #expand("02_hybrid_assembly/{assembly_sample}/{assembly_sample}_assembly.fasta", assembly_sample=SAMP),
        expand("03_indexed_anno_reference/{assembly_sample}/{assembly_sample}.gbk",  assembly_sample=SAMP),
        expand("04_mlst/{assembly_sample}_mlst", assembly_sample=SAMP),
        expand("05_assembly_blast_databases/{assembly_sample}_assembly.fasta", assembly_sample=SAMP),
        expand("05_assembly_blast_databases/{assembly_sample}_assembly.fasta.{extension}", assembly_sample=SAMP, extension=["nhr","nin","nog","nsd","nsi","nsq"]),
        expand("06_blast_all_bacteriocins/{assembly_sample}_blastoutput_all_bacteriocins.txt", assembly_sample=SAMP),
        expand("07_blast_enterocinA/{assembly_sample}_blastoutput_enterocinA.txt", assembly_sample=SAMP),
        expand("08_blast_bac43/{assembly_sample}_blastoutput_bac43.txt", assembly_sample=SAMP),
        expand("09_blast_ent53B/{assembly_sample}_blastoutput_ent53B.txt", assembly_sample=SAMP),
        expand("10_blast_vanABMD/{assembly_sample}_blastoutput_vanABMD.txt", assembly_sample=SAMP),
        expand("11_snippy/{reference}/{sample}/snps.consensus.fa", sample=SAMP, reference=REFERENCE),
        expand("12_core_genome/{reference}.full.aln", reference=REFERENCE),
        expand("12_core_genome/{reference}.clean.full.aln", reference=REFERENCE),
        expand("13_gubbins/{reference}.final_tree.tre", reference=REFERENCE),
        expand("13_gubbins/{reference}.node_labelled.final_tree.tre", reference=REFERENCE),



#Trimmoamatic
#File should be format ######_R1_001.fastq.gz 
rule trimmomatic_pe:
    input:
        r1="shortread/{sample}_R1_001.fastq.gz",
        r2="shortread/{sample}_R2_001.fastq.gz"
    output:
        r1="01_trimmed_shortread/{sample}_01.fq.gz",
        r2="01_trimmed_shortread/{sample}_02.fq.gz",
        r1_unpaired="01_trimmed_shortread/{sample}_01_unpaired.fq.gz",
        r2_unpaired="01_trimmed_shortread/{sample}_02_unpaired.fq.gz"
    conda:
        "envs/trimmomatic_2.yml"
    params:
        bp=config["Trimmomatic"]
    threads: 4
    log: "logs/{sample}_trimmomatic.log"
    shell:"""
    touch logs
    trimmomatic PE -phred33 -threads {threads} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} {params.bp}    2> {log}
    fastqc {output.r1} -o logs
    fastqc {output.r2} -o logs
    """

#trim long reads
rule demultiplex_filter_nanopore:
    input:
        r3="longread/{longread_sample}_merged_LR.fastq"
    output:
        r4="01_trimmed_longread/{longread_sample}_longread.fq",
        temp=temp("01_trimmed_longread/{longread_sample}_temp.fq")
    log: log1="logs/{longread_sample}_nanofilt.log", log2="logs/{longread_sample}_porechop.log", log3="logs/{longread_sample}_nanoplot"
    params:
        quality=config["Nanofilt"]["quality"],
        length=config["Nanofilt"]["length"],
        adapter_threshold=config["Porechop"]["adapter_threshold"]
    threads: 12
    conda:
        "envs/long_read_qc_3.yml"
    shell:"""
    touch logs
    porechop -i {input.r3} -o {output.temp} -t {threads} --adapter_threshold {params.adapter_threshold} 2> {log.log2}
    cat {output.temp} | NanoFilt -q {params.quality} -l {params.length} --logfile {log.log1} > {output.r4}
    #NanoPlot --fastq {output.r4} --N50 -o {log.log3}
    """
    
# #de novo assembly 
rule short_read_assembly:
    input:
        short_read_01="01_trimmed_shortread/{assembly_sample}_01.fq.gz",
        short_read_02="01_trimmed_shortread/{assembly_sample}_02.fq.gz",
    output:
        assembly_file="02_assembly/{assembly_sample}/{assembly_sample}_assembly.fasta",
        log_file="02_assembly/{assembly_sample}/{assembly_sample}.log"
    params:
        assembly_dir="02_assembly/{assembly_sample}",
        mode=config["Unicycler"]["mode"]
    threads: 32
    log: "logs/{assembly_sample}_unicycler.log"
    conda: "envs/unicycler_5.yml"
    priority: 100
    shell:"""
         unicycler -1 {input.short_read_01} -2 {input.short_read_02} -o {params.assembly_dir} -t {threads} --mode {params.mode}
         mv {params.assembly_dir}/unicycler.log {output.log_file}
         mv {params.assembly_dir}/assembly.fasta {output.assembly_file}
         cp {output.log_file} {log}
    """

#Specify hybrid assembly fasta file as output
rule long_read_assembly:
    input:
        short_read_01="01_trimmed_shortread/{assembly_sample}_01.fq.gz",
        short_read_02="01_trimmed_shortread/{assembly_sample}_02.fq.gz",
        long_reads="01_trimmed_longread/{assembly_sample}_longread.fq"
    output:
        assembly_file="02_hybrid_assembly/{assembly_sample}/{assembly_sample}_assembly.fasta",
        log_file="02_hybrid_assembly/{assembly_sample}/{assembly_sample}.log"
    params:
        assembly_dir="02_hybrid_assembly/{assembly_sample}",
        mode=config["Unicycler"]["mode"]
    threads: 8
    log: "logs/{assembly_sample}_unicycler.log"
    conda: "envs/unicycler_4.yml"
    priority: 100
    shell:"""
         unicycler -1 {input.short_read_01} -2 {input.short_read_02} -l {input.long_reads} -o {params.assembly_dir} -t {threads} --mode {params.mode}
         mv {params.assembly_dir}/unicycler.log {output.log_file}
         mv {params.assembly_dir}/assembly.fasta {output.assembly_file}
         cp {output.log_file} {log}
    """

rule index_reference:
    input:
        reference_fasta="02_assembly/{assembly_sample}/{assembly_sample}_assembly.fasta"
    output:
        reference_fasta="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}_assembly.fasta",
        indexed_reference=expand("03_indexed_anno_reference/{{assembly_sample}}/{{assembly_sample}}_assembly.fasta.{extension}", extension=INDEX_FILE),
        dictionary="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}_assembly.dict"
    params:
        output_dir="03_indexed_anno_reference/{assembly_sample}",
        input_file="{assembly_sample}_assembly.fasta"
    log: "logs/{assembly_sample}_unicycler.log"
    threads: 12
    conda:
         "envs/index_mapping_1.yml"
    shell:"""
    cp {input.reference_fasta} {params.output_dir}
    cd {params.output_dir}
    bwa index -a bwtsw {params.input_file}
    samtools faidx {params.input_file}
    picard CreateSequenceDictionary REFERENCE={params.input_file}
    bwa index {params.input_file}
    """

rule annotate_reference:
    input:
        reference_fasta="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}_assembly.fasta"
    output:
        gff_annotation="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}.gff",
        gbf_annotation="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}.gbk"
    conda:
        "envs/prokka_1.yaml"
    params:
        prefix="{assembly_sample}",
        gff_loc="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}/{assembly_sample}.gff",
        gbf_loc="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}/{assembly_sample}.gbk",
        prokka_dir="03_indexed_anno_reference/{assembly_sample}/{assembly_sample}",
        prokka_dir_2="03_indexed_anno_reference/{assembly_sample}"
    threads: 12
    shell:"""
    prokka {input.reference_fasta} --prefix {params.prefix} --outdir {params.prokka_dir} --force
    cp {params.gff_loc} {params.prokka_dir_2}
    cp {params.gbf_loc} {params.prokka_dir_2}
    rm -r {params.prokka_dir}
    """

rule identify_mlst:
    input:
        assembly_file="02_assembly/{assembly_sample}/{assembly_sample}_assembly.fasta"
    output:
        mlst_file="04_mlst/{assembly_sample}_mlst",
    conda: "envs/mlst_1.yml"
    params: 
        output_dir = "04_mlst"
    threads: 6
    shell:"""
        mlst {input.assembly_file} > {output.mlst_file}
    """


rule make_blast_databases:
    input:
        assembly_file = "02_assembly/{assembly_sample}/{assembly_sample}_assembly.fasta"
    output:
        assembly = "05_assembly_blast_databases/{assembly_sample}_assembly.fasta",
        blast_db = multiext("05_assembly_blast_databases/{assembly_sample}_assembly.fasta", ".nhr",".nin",".nog",".nsd",".nsi",".nsq")
    conda: "envs/blast_1.yml"
    params:
        output_dir = "05_assembly_blast_databases",
    threads: 12
    shell: """
        cp {input.assembly_file} {params.output_dir}
        makeblastdb -in {output.assembly} -input_type fasta -parse_seqids -dbtype nucl
    """

rule blast_bacteriocins:
    input:
        query = "blast_fastas/all_bacteriocin_classes.fa",
        database = "05_assembly_blast_databases/{assembly_sample}_assembly.fasta",
        blast_files = multiext("05_assembly_blast_databases/{assembly_sample}_assembly.fasta", ".nhr",".nin",".nog",".nsd",".nsi",".nsq")
    output:
        bacteriocin_hits = "06_blast_all_bacteriocins/{assembly_sample}_blastoutput_all_bacteriocins.txt"
    conda: "envs/blast_1.yml"
    threads: 12
    shell: """
        tblastn -db {input.database} -query {input.query} -out {output.bacteriocin_hits} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" 
    """

rule blast_enterocinA:
    input:
        query = "blast_fastas/enterocinA.fa",
        database = "05_assembly_blast_databases/{assembly_sample}_assembly.fasta"
    output:
        enterocinA_hits = "07_blast_enterocinA/{assembly_sample}_blastoutput_enterocinA.txt"
    conda: "envs/blast_1.yml"
    threads: 8
    shell: """
        blastn -db {input.database} -query {input.query} -out {output.enterocinA_hits} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" 
    """

rule blast_bac43:
    input:
        query = "blast_fastas/bac43.fa",
        database = "05_assembly_blast_databases/{assembly_sample}_assembly.fasta"
    output:
        bac43_hits = "08_blast_bac43/{assembly_sample}_blastoutput_bac43.txt"
    conda: "envs/blast_1.yml"
    threads: 8
    shell: """
        blastn -db {input.database} -query {input.query} -out {output.bac43_hits} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" 
    """

rule blast_ent53B:
    input:
        query = "blast_fastas/ent53B.fa",
        database = "05_assembly_blast_databases/{assembly_sample}_assembly.fasta"
    output:
        ent53B_hits = "09_blast_ent53B/{assembly_sample}_blastoutput_ent53B.txt"
    conda: "envs/blast_1.yml"
    threads: 8
    shell: """
        blastn -db {input.database} -query {input.query} -out {output.ent53B_hits} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" 
    """

rule blast_vanABMD:
    input:
        query = "blast_fastas/vanABMD.fa",
        database = "05_assembly_blast_databases/{assembly_sample}_assembly.fasta"
    output:
        van_hits = "10_blast_vanABMD/{assembly_sample}_blastoutput_vanABMD.txt"
    conda: "envs/blast_1.yml"
    threads: 8
    shell: """
        blastn -db {input.database} -query {input.query} -out {output.van_hits} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" 
    """

rule snippy:
    input:
        r1="01_trimmed_shortread/{sample}_01.fq.gz",
        r2="01_trimmed_shortread/{sample}_02.fq.gz",
        reference="snippy_reference/{reference}.gbk"
    output:
        consensus="11_snippy/{reference}/{sample}/snps.consensus.fa",
        consensus_dir=directory("11_snippy/{reference}/{sample}/")
    params:
        id="{sample}",
        dir="11_snippy/{reference}/{sample}/"
    conda:
        "envs/snippy_4.yaml"
    threads: 8
    shell:"""
    snippy --outdir {params.dir} --R1 {input.r1} --R2 {input.r2} --ref {input.reference} --cpus {threads} --mincov 20 --minfrac .9 --mapqual 60 --rgid {params.id} --force
    """

rule snippy_core:
    input:
        dirs=expand("11_snippy/{reference}/{sample}", reference=REFERENCE, sample=SAMP_CORE),
        reference="snippy_reference/{reference}.gbk"
    output:
        aln="12_core_genome/{reference}.full.aln",
        clean_aln="12_core_genome/{reference}.clean.full.aln"
    params:
        id="{reference}"
    conda:
        "envs/snippy_4.yaml"
    shell:"""
    snippy-core --prefix=12_core_genome/{params.id} --ref={input.reference} {input.dirs}
    snippy-clean_full_aln 12_core_genome/{params.id}.full.aln > 12_core_genome/{params.id}.clean.full.aln
    """

rule gubbins:
    input:
        clean_aln="12_core_genome/{reference}.clean.full.aln"
    output:
        final_tree = "13_gubbins/{reference}.final_tree.tre",
        labeled_final_tree = "13_gubbins/{reference}.node_labelled.final_tree.tre"
    params:
        prefix = "13_gubbins/{reference}"
    conda: "envs/gubbins.yml"
    threads: 24
    shell: """
    run_gubbins.py -v -p {params.prefix} --threads {threads} {input.clean_aln} --filter-percentage 25 --first-tree-builder rapidnj --first-model JC --tree-builder raxmlng --model GTR 
    """


