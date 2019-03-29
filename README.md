## Introduction
------------------------------
This is an attempt to produce the nullarbor report (https://github.com/tseemann/nullarbor) using snakemake (https://snakemake.readthedocs.io/en/stable/), hence the name `snakeNullarbor`   

### To-Do
* allow each of the snake files to run independently if needed   

### Current limitation
* Visualising the phylogentic tree and plotting metdata is done outside Snakemake  
* The option `--use-conda` within snakemake is not currently feasible 
* The config file **MUST** include both fastq and fasta files 

### Dependencies 
------------------------------
The workflow was used with the following versions of software   
* SPAdes v3.12.0 http://cab.spbu.ru/software/spades/ 
* Kraken v0.10.6 https://ccb.jhu.edu/software/kraken/
* Prokka v1.13.3 https://github.com/tseemann/prokka
* Roary v3.6.1 https://sanger-pathogens.github.io/Roary/
* Bowtie2 v2.3.0 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* BWA v0.7.12-r1039 http://bio-bwa.sourceforge.net/
* SAMTools v1.4-14-g90b995f  http://samtools.sourceforge.net/
* FastTree v2.1.9 http://www.microbesonline.org/fasttree/

## Getting Started
------------------------------
### Obtaining the Current Version of the workflow 
* Set up a project folder for the run 

* Download the latest version from github  
`git clone https://gitlab.com/Mostafa.Abdel-Glil/snakenullarbor.git`

### generate a config file for the pipeline 

### Running the snake pipeline 
 
`snakemake -np --quiet --snakefile master.Snakefile`

## Run a test dataset

## Contact 
___
comments should be adderessed to Mostafa.Abdel-Glil@fli.de 
