"""----------------------------
Friedrich-Loeffler-Institut (https://www.fli.de/), IBIZ
date: March, 20, 2019
Author: Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)
-------------------------------
This is an attempt to produce the nullarbor report (https://github.com/tseemann/nullarbor) using snakemake (https://snakemake.readthedocs.io/en/stable/), hence the name `snakeNullarbor`
-------------------------------

# TODO:
* isolates.txt, use python code to prepare it nicely
* make the log files more nice and clear
* in the MLST and snp-dists rule, results folder in the sed command is hardcoded

Assembly
* assembler_options =config["assembler_options"] ---snake variable and rule
* -opts {params.assembler} options does not wotŕk in the bash --- bash
 * shovill --outdir $OUTDIR --cpus $THREADS --R1 $READ1 --R2 $READ2 --force --ram 10 $assembler_options $assembler_options_additional

Prokka
* --gram pos can not be used, no `signalp` in conda

Nullarbor
* remove the tools versions in nullarbor and make it read from the envs folder
* export KRAKEN_DEFAULT_DB=/home/mostafa.abdel/dbs/miniKraken/minikraken_20171013_4GB for DBs, or better read from the config file
"""
import os
import tempfile
configfile: "config.yaml"
#include: snakefiles + "defineFolders.Snakefile"
#defineFolders
#working_dir=config['working_dir']
snakemake_folder=config['snakemake_folder']
raw_data_dir=config['raw_data_dir']
fasta_dir=config['fasta_dir']
results_dir="results/" #hardcoded in certain rules, don't change
temporary_todelete= results_dir + "temporary_todelete/"

#samples
SAMPLES, = glob_wildcards( raw_data_dir + "{sample}_R1.fastq.gz")
GENOMES, = glob_wildcards( fasta_dir + "{genome}.fasta")
SAMPLES = set(SAMPLES)
GENOMES = set(GENOMES)
DATASET =  [SAMPLES, GENOMES]
#scripts_paths
bin_dir= snakemake_folder + "bin/",
lib_dir= snakemake_folder + "lib/"
nullarbor_bin= snakemake_folder + "lib/nullarbor/bin/"
#envs_folder=config["bin_dir"] + "../lib/envs/"
envs_folder="lib/envs/"
#variables
if "AMR_db_abricate" in config:
    AMR_db_abricate= config["AMR_db_abricate"]
else:
    AMR_db_abricate="ncbi"
if "VF_db_abricate" in config:
    VF_db_abricate= config["VF_db_abricate"]
else:
    VF_db_abricate="vfdb"
if "prokka_genustag" in config:
    prokka_genustag= config["prokka_genustag"]
else:
    prokka_genustag= ""
if "prokka_params" in config:
    prokka_params= config["prokka_params"]
else:
    prokka_params= ""
#mlst_parameters
if "MLST_options" in config:
    MLST_options= config["MLST_options"]
else:
    MLST_options= ""
assembler= "shovill",
assembler_options= ""
taxoner= config["taxoner"]
reference= config["reference"]

rule all:
    input:
        seq_yield = expand(  results_dir + "{sample}/yield.tab", sample=SAMPLES),
        Assembly =expand( results_dir + "{sample}/contigs.fa", sample=SAMPLES),
        kraken_tab= expand(results_dir + "{sample}/kraken.tab", sample=SAMPLES),
        kraken_tab_assemblies= expand(results_dir + "{genome}/kraken.tab", genome=GENOMES),
        virulome =expand( results_dir + "{sample}/virulome.tab", sample=SAMPLES),
        virulome_assemblies =expand( results_dir + "{genome}/virulome.tab", genome=GENOMES),
        resistome =expand( results_dir + "{sample}/resistome.tab", sample=SAMPLES),
        resistome_assemblies =expand( results_dir + "{genome}/resistome.tab", genome=GENOMES),
        contig = expand(  results_dir + "{genome}/contigs.fa", genome=GENOMES),
        #prokka_gff= expand(results_dir + "{sample}/prokka/{sample}.gff", sample=SAMPLES),
        #prokka_assemblies_gff= expand( results_dir + "{genome}/prokka/{genome}.gff", genome=GENOMES),
        isolates_list= results_dir + "isolates.txt",
        denovo = results_dir + "denovo.tab",
        snippy_snps= expand( results_dir + "{sample}/snippy/snps.tab", sample=SAMPLES),
        snippy_snps_assemblies= expand( results_dir + "{genome}/snippy/snps.tab", genome=GENOMES),
        snippycore= results_dir + "core.aln",
        Roary= results_dir + "roary/pan_genome_reference.fa",
        Roary_dir=results_dir + "roary/pangenome_frequency.png",
        phylogeny_tree= results_dir + "core.newick",
        pan_svg= results_dir + "roary/pan.svg",
        mlstresults= results_dir + "mlst.tab",
        snpdistsresults= results_dir + "distances.tab",
        nullarborreport= results_dir + "report/index.html",

"""
collect assemblies
"""
rule collect_assemblies:
    input:
        assemblies = fasta_dir + "{genome}.fasta",
        yield_na = lib_dir + "yield.tab"
    output:
        contig = results_dir + "{genome}/contigs.fa",
        contig_yield = results_dir + "{genome}/yield.tab",
    conda:
        envs_folder + "bioawk.yaml" #spades, sickle
    shell:
        "bash {bin_dir}fa_rename.sh {input.assemblies} {output.contig} && cp -v -f {input.yield_na} {output.contig_yield}" #ln -s -f ##* Prokka does not like contigs ID > 37, '--centre X --compliant' is not the way always
"""
assembly
"""
rule assembly:
    input:
        r1 = raw_data_dir + "{sample}_R1.fastq.gz",
        r2 = raw_data_dir + "{sample}_R2.fastq.gz",
    output:
        contig = results_dir + "{sample}/contigs.fa",
    threads: 16 #increasing threads produces errors with spades
    conda:
        envs_folder + "shovill.yaml" #spades, sickle
    params:
        assembler = assembler_options,
        assembly_dir = directory (results_dir + "{sample}")
    shell:
        "bash {bin_dir}assembly.sh -a {assembler} --R1 {input.r1} --R2 {input.r2} -o {output.contig} -d {params.assembly_dir}/{assembler}  -p {threads} -f true " #-t trimmomatic
"""
kraken
"""
rule kraken:
    input:
        r1 = raw_data_dir + "{sample}_R1.fastq.gz",
        r2 = raw_data_dir + "{sample}_R2.fastq.gz",
        db=config["kraken"],
        #taxoner= config["taxoner"]
    output:
        kraken_tab= results_dir + "{sample}/kraken.tab"
    threads: 32
    conda:
      envs_folder + "kraken.yaml"
    shell:
        "kraken --db {input.db} --paired --check-names --threads {threads} --gzip-compressed --fastq-input {input.r1} {input.r2} | kraken-report --db {input.db} > {output.kraken_tab} 2>&1 | sed 's/^/[kraken] /'"
rule kraken_assemblies:
    input:
        assemblies = fasta_dir + "{genome}.fasta",
        db=config["kraken"],
        #taxoner= config["taxoner"]
    output:
        kraken_tab= results_dir + "{genome}/kraken.tab"
    threads: 32
    conda:
      envs_folder + "kraken.yaml"
    shell:
        "kraken --db {input.db} --threads {threads} --fasta-input {input.assemblies} | kraken-report --db {input.db} > {output.kraken_tab} 2>&1 | sed 's/^/[kraken] /'"
"""
denovo fasta
"""
rule fasta_denovo:
    input:
        fasta = expand( results_dir + "{sample}/contigs.fa", sample=SAMPLES),
        fasta_assembly = expand( results_dir + "{sample}/contigs.fa", sample=GENOMES),
    output:
        denovo = results_dir + "denovo.tab",
    conda:
        envs_folder + "perl.yaml"
    shell:
        "{nullarbor_bin}fa --minsize 500 -e -t {input.fasta} {input.fasta_assembly} | sed 's#results/##g' | tee -a {output.denovo}"
"""
reference
"""
rule reference:
    output:
        ref_fasta= results_dir + "ref.fa",
    conda:
        envs_folder + "reference.yaml"
    shell:
        "seqret -auto -filter -osformat2 fasta {reference} > {output.ref_fasta} && \
        samtools faidx {output.ref_fasta}"
"""
isolates.txt #revision
"""
names =expand( "{sample}", sample=SAMPLES),
names_genomes =expand( "{sample}", sample=GENOMES),
isolates = [names, names_genomes]
rule isolates_list:
    output:
        isolates_list= results_dir + "isolates.txt",
        tempfile=temp(results_dir + "tempfile",)
    shell:
        '''
        echo -ne "{isolates}" > {output.tempfile} && bash {bin_dir}fixnames.sh {output.tempfile} > {output.isolates_list}
        '''
#import sys
#sys.stdout = open('results/isolates.txt', 'w')
#for w in isolates:
#print(w, sep="\n")
#print (*names, sep='\n')
#print (*names_genomes, sep='\n')
"""
yield sequences
"""
rule seq_yield:
    input:
        r1 = raw_data_dir + "{sample}_R1.fastq.gz",
        r2 = raw_data_dir + "{sample}_R2.fastq.gz",
        ref_fasta= results_dir + "ref.fa",
    output:
        seq_yield = results_dir + "{sample}/yield.tab",
    conda:
        envs_folder + "perl.yaml"
    shell:
        "{nullarbor_bin}fq --quiet --ref {input.ref_fasta} {input.r1} {input.r2} > {output.seq_yield}"
"""
virulome
"""
rule virulome:
    input:
        contig = results_dir + "{sample}/contigs.fa",
    output:
        virulome = results_dir + "{sample}/virulome.tab",
    threads: 16
    conda:
        envs_folder + "abricate.yaml"
    params:
        VF_db_abricate = VF_db_abricate,
    shell:
        "abricate --db {params.VF_db_abricate} --threads {threads} {input.contig} > {output.virulome}"

"""
resistome
"""
rule resistome:
    input:
        contig = results_dir + "{sample}/contigs.fa",
    output:
        resistome = results_dir + "{sample}/resistome.tab",
    threads: 16
    conda:
        envs_folder + "abricate.yaml"
    params:
        AMR_db_abricate = AMR_db_abricate,
    shell:
        "abricate --db {params.AMR_db_abricate} --threads {threads} {input.contig} > {output.resistome}"

"""
Prokka
"""
rule prokka:
    input:
        contig= results_dir + "{sample}/contigs.fa",
    output:
        prokka_gff= results_dir + "{sample}/prokka/{sample}.gff",
        prokka_gbk= results_dir + "{sample}/prokka/{sample}.gbk",
    threads: 32
    conda: envs_folder + "prokka.yaml"
    params:
        prokka_outdir= results_dir + "{sample}/prokka",
        prefix="{sample}",
        strain_name="{sample}",
        locustag= "{sample}",
    shell:
        "prokka --cpus {threads} --force --kingdom Bacteria --outdir {params.prokka_outdir}/ --prefix {params.prefix} --genus {prokka_genustag} --gcode 11 --locustag {params.locustag} --strain {params.strain_name} {prokka_params} {input.contig} 2>&1 | sed 's/^/[prokka] /' " #--centre Institute --compliant to generate clean contig names
rule prokka_ln:
    input:
        prokka_gff= results_dir + "{sample}/prokka/{sample}.gff",
        prokka_gbk= results_dir + "{sample}/prokka/{sample}.gbk",
    output:
        prokka_gff= results_dir + "{sample}/contigs.gff",
        prokka_gbk= results_dir + "{sample}/contigs.gbk",
    params:
        prokka_outdir= results_dir + "{sample}/prokka",
    shell:
        "ln -s -f $(realpath {input.prokka_gff}) {output.prokka_gff}; ln -s -f $(realpath {input.prokka_gbk}) {output.prokka_gbk}"
"""
snippy
"""
#export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64
rule snippy:
    input:
        r1 = raw_data_dir + "{sample}_R1.fastq.gz",
        r2 = raw_data_dir + "{sample}_R2.fastq.gz",
    output:
        snippy_snps= results_dir + "{sample}/snippy/snps.tab",
        #snippy_folder=  results_dir + "{sample}",
    threads: 32
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_dir + "{sample}/snippy",
        #snippy_outdir_cp= results_dir + "{sample}/",
    shell:
        "snippy --force  --cpus {threads} --ram {threads} --outdir {params.snippy_outdir} --ref {reference} --R1 {input.r1} --R2 {input.r2} 2>&1 | sed 's/^/[snippy] /'"
rule snippy_ln:
    input:
        snippy_snps= results_dir + "{sample}/snippy/snps.tab",
    output:
        snippy_snps_tab=  results_dir + "{sample}/snps.tab",
        snippy_snps_aligned= results_dir + "{sample}/snps.aligned.fa",
        snippy_snps_rawvcf= results_dir + "{sample}/snps.raw.vcf",
        snippy_snps_vcf= results_dir + "{sample}/snps.vcf",
        snippy_snps_bam= results_dir + "{sample}/snps.bam",
        snippy_snps_bai= results_dir + "{sample}/snps.bam.bai",
        snippy_snps_log= results_dir + "{sample}/snps.log",
#    threads: 32
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_dir + "{sample}/snippy",
    shell:
        "ln -s -f $(realpath {params.snippy_outdir}/snps.tab) {output.snippy_snps_tab}; ln -s -f $(realpath {params.snippy_outdir}/snps.aligned.fa) {output.snippy_snps_aligned} ; ln -s -f $(realpath {params.snippy_outdir}/snps.raw.vcf) {output.snippy_snps_rawvcf}; ln -s -f $(realpath {params.snippy_outdir}/snps.vcf) {output.snippy_snps_vcf}; ln -s -f $(realpath {params.snippy_outdir}/snps.bam) {output.snippy_snps_bam}; ln -s -f $(realpath {params.snippy_outdir}/snps.bam.bai) {output.snippy_snps_bai}; ln -s -f $(realpath {params.snippy_outdir}/snps.log) {output.snippy_snps_log}"
rule snippy_assemblies:
    input:
        contig = fasta_dir + "{genome}.fasta"
    output:
        snippy_snps= results_dir + "{genome}/snippy/snps.tab",
    threads: 32
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_dir + "{genome}/snippy",
    shell:
        "snippy --force  --cpus {threads} --ram {threads} --outdir {params.snippy_outdir} --ref {reference} --ctgs {input.contig}  2>&1 | sed 's/^/[snippy] /'"

def snippy_folders(wildcards):
    return expand(results_dir + "{sample}", sample=SAMPLES)
def snippy_folders_assemblies(wildcards):
    return expand(results_dir + "{genome}", genome=GENOMES)

rule snippy_core:
    input:
        #snippy_folders=  expand( results_dir + "{sample}", sample=SAMPLES),
        #snippy_folders_assemblies=  expand( results_dir + "{genome}", genome=GENOMES),
        a=expand( results_dir + "{sample}/snps.tab", sample=SAMPLES),
        b=expand( results_dir + "{sample}/snps.aligned.fa", sample=SAMPLES),
        c=expand( results_dir + "{sample}/snps.raw.vcf", sample=SAMPLES),
        d=expand( results_dir + "{sample}/snps.vcf", sample=SAMPLES),
        e=expand( results_dir + "{sample}/snps.bam", sample=SAMPLES),
        f=expand( results_dir + "{sample}/snps.bam.bai", sample=SAMPLES),
        g=expand( results_dir + "{sample}/snps.log", sample=SAMPLES),
        h=expand( results_dir + "{genome}/snps.tab", genome=GENOMES),
        i=expand( results_dir + "{genome}/snps.aligned.fa", genome=GENOMES),
        j=expand( results_dir + "{genome}/snps.raw.vcf", genome=GENOMES),
        k=expand( results_dir + "{genome}/snps.vcf", genome=GENOMES),
        l=expand( results_dir + "{genome}/snps.bam", genome=GENOMES),
        m=expand( results_dir + "{genome}/snps.bam.bai", genome=GENOMES),
        n=expand( results_dir + "{genome}/snps.log", genome=GENOMES),
    output:
        snippycore= results_dir + "core.aln",
    conda: envs_folder + "snippy.yaml"
    params:
        snippycoreoutdir= results_dir + "core",
        snippy_folders = snippy_folders,
        snippy_folders_assemblies = snippy_folders_assemblies,
        #snippy_outdir_cp= results_dir + "{sample}/",
    shell:
        "snippy-core  --ref {reference} {params.snippy_folders} {params.snippy_folders_assemblies} --prefix {params.snippycoreoutdir}  2>&1 | sed 's/^/[snippy-core] /'"
"""
Core genome phylogeny using FastTree
"""
rule FastTree: #run fasttree
    input:
        snippycore= results_dir + "core.aln",
    output:
        phylogeny_tree= results_dir + "core.newick",
        phylogeny_svg= results_dir + "core.svg"
    threads: 64
    #benchmark: benchmarks_folder + "FastTree.txt"
    conda: envs_folder + "FastTree.yaml" #needs revision
    shell:
        "FastTree -nt -gtr {input.snippycore} > {output.phylogeny_tree} && \
        nw_display -S -s -w 1024 -l 'font-size:12;font-family:sans-serif;' -i 'opacity:0' -b 'opacity:0' -v 16 {output.phylogeny_tree}  > {output.phylogeny_svg}"

"""
Roary
"""
#export PERL5LIB=$PERL5LIB:/home/software/Roary/lib/
rule Roary: #run roary
    input:
        #bin_dir=config["bin_dir"],
        prokka_gff= expand(results_dir + "{sample}/prokka/{sample}.gff", sample=SAMPLES),
        prokka_assemblies_gff= expand(results_dir + "{genome}/prokka/{genome}.gff", genome=GENOMES),
    output:
        tmp_dir=temp(directory(temporary_todelete)),
        Roary_pangenome_fa= results_dir + "roary/pan_genome_reference.fa",
        Roary_aln= results_dir + "roary/core_gene_alignment.aln",
        roary_presenceabsence= results_dir + "roary/gene_presence_absence.csv",
        roary_acc= results_dir + "roary/accessory_binary_genes.fa.newick",
    threads: 64
    log: results_dir + "roary.log"
    conda: envs_folder + "roary.yaml"
    params:
        options=config["roary_params"],
        Roary_dir=results_dir + "roary"
    shell:
        "bash {bin_dir}fixRoaryOutDirError.sh {params.Roary_dir} {output.tmp_dir}"
        " && roary -p {threads} -f {params.Roary_dir} {params.options} {input.prokka_gff} {input.prokka_assemblies_gff} 2>&1 | sed 's/^/[roary] /' | tee -a {log}" #-e, core genes alignment using PRANK, -r, Rplots ,
        # sometimes you may need to set -i (minimum percentage identity for blastp) and -s (dont split paralogs), according to the organism
        #" && python3 {bin_dir}roary_plots.py {params.Roary_dir}/accessory_binary_genes.fa.newick {params.Roary_dir}/gene_presence_absence.csv"
        #" && mv -t {params.Roary_dir} pangenome_frequency.png pangenome_matrix.png pangenome_pie.png"
rule Roary_plots:
    input:
        roary_presenceabsence= results_dir + "roary/gene_presence_absence.csv",
        roary_acc= results_dir + "roary/accessory_binary_genes.fa.newick",
    output:
        Roary_dir=results_dir + "roary/pangenome_frequency.png"
    conda: envs_folder + "roary_plots.yaml"
    params:
        #options=config["roary_params"],
        Roary_dir=results_dir + "roary"
    shell:
        " python3 {bin_dir}roary_plots.py {params.Roary_dir}/accessory_binary_genes.fa.newick {params.Roary_dir}/gene_presence_absence.csv"
        " && mv -t {params.Roary_dir} pangenome_frequency.png pangenome_matrix.png pangenome_pie.png"
#source ~/.bash_profile SVG.pm
rule Roary_svg: #run roary
    input:
        roary_pangenome_fa= results_dir + "roary/pan_genome_reference.fa",
        roary_presenceabsence= results_dir + "roary/gene_presence_absence.csv",
        roary_acc= results_dir + "roary/accessory_binary_genes.fa.newick",
    output:
        acc_svg= results_dir + "roary/acc.svg",
        pan_svg= results_dir + "roary/pan.svg",
    #threads: 64
    conda: envs_folder + "perl.yaml"
    params:
        options=config["roary_params"],
        Roary_dir=results_dir + "Roary"
    shell:
        " nw_display -S -s -w 1024 -l 'font-size:12;font-family:sans-serif;' -i 'opacity:0' -b 'opacity:0' -v 16 {input.roary_acc} > {output.acc_svg}"
        " && {nullarbor_bin}roary2svg.pl {input.roary_presenceabsence} > {output.pan_svg}"

"""
mlst
"""
rule MLST:
    input:
        contigs= expand(results_dir + "{sample}/contigs.fa", sample=SAMPLES),
        contigs_assemblies= expand(results_dir + "{genome}/contigs.fa", genome=GENOMES),
        ref_fasta= results_dir + "ref.fa",
    output:
        mlstresults= results_dir + "mlst.tab",
    conda:
        envs_folder + "mlst.yaml"
    shell:
        "mlst --quiet {MLST_options} {input.ref_fasta} {input.contigs} {input.contigs_assemblies} | sed 's#results/##g' | tee -a {output.mlstresults}"
"""
snpdists
"""
rule snpdists:
    input:
        snippycore= results_dir + "core.aln",
    output:
        snpdistsresults= results_dir + "distances.tab",
    conda:
        envs_folder + "snpdists.yaml"
    shell:
        "snp-dists -b {input.snippycore} > {output.snpdistsresults}"
"""
Nullarbor report
"""
rule Nullarbor:
    input:
        seq_yield = expand(  results_dir + "{sample}/yield.tab", sample=SAMPLES),
        Assembly =expand( results_dir + "{sample}/contigs.fa", sample=SAMPLES),
        virulome =expand( results_dir + "{sample}/virulome.tab", sample=SAMPLES),
        virulome_assemblies =expand( results_dir + "{genome}/virulome.tab", genome=GENOMES),
        resistome =expand( results_dir + "{sample}/resistome.tab", sample=SAMPLES),
        resistome_assemblies =expand( results_dir + "{genome}/resistome.tab", genome=GENOMES),
        contig = expand(  results_dir + "{genome}/contigs.fa", genome=GENOMES),
        prokka_gff= expand(results_dir + "{sample}/contigs.gff", sample=SAMPLES),
        prokka_assemblies_gff= expand(results_dir + "{genome}/contigs.gff",  genome=GENOMES),
        #prokka_assemblies_gff= expand( results_dir + "{genome}/prokka/{genome}.gff", genome=GENOMES),
        isolates_list= results_dir + "isolates.txt",
        denovo = results_dir + "denovo.tab",
        #snippy_snps= expand( results_dir + "{sample}/snippy/snps.tab", sample=SAMPLES),
        #snippy_snps_assemblies= expand( results_dir + "{genome}/snippy/snps.tab", genome=GENOMES),
        snippycore= results_dir + "core.aln",
        Roary= results_dir + "roary/pan_genome_reference.fa",
        phylogeny_tree= results_dir + "core.newick",
        pan_svg= results_dir + "roary/pan.svg",
        mlstresults= results_dir + "mlst.tab",
        snpdistsresults= results_dir + "distances.tab",
    output:
        nullarborreport= results_dir + "report/index.html",
    params:
        nullarborfolder = "report",
    conda:
        envs_folder + "perl.yaml"
    shell:
        "cd {results_dir} && {nullarbor_bin}nullarbor-report.pl --name SampleReport --indir . --outdir {params.nullarborfolder}"
        " && mkdir -p $HOME/public_html/MDU/PROJNAME"
        " && install -p -D -t $HOME/public_html/MDU/PROJNAME report/*"
