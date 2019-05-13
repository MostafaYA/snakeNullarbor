## Introduction
------------------------------
This is an attempt to produce the [nullarbor report](https://github.com/tseemann/nullarbor) using [snakemake](https://snakemake.readthedocs.io/en/stable/)

<img src="workflowpic.png" width="1000" />

## Getting Started

### download the project 
* Set up a project folder for the run 
cd /path/to/installation
* Download the latest version from gitlab  
`git clone https://gitlab.com/FLI_Bioinfo/snakenullarbor.git`

## Running snakenullarbor 

### 1. Run the bash script `snakenullarbor.sh`
To start snakenullarbor, run the bash script `snakenullarbor.sh`. Detailed usage is below 

##### Dependencies 
* pigz (http://zlib.net/pigz/)
* seqret (http://emboss.sourceforge.net/apps/cvs/emboss/apps/seqret.html)

```
USAGE:
   bash ./snakenullarbor.sh -fq fastq_directory -r REF
REQUIRED:
   -fq, --fastq-directory DIR
   -r, --ref 'Reference_strain'
OPTIONAL:
   -fa, --fasta-directory DIR
   -o, --outdir (default: input/)
   -h, --help

```


* As input for `snakenullarbor.sh`, you need to use the folder where the fastq files are present with option `-fq`. In standard cases, fastq files could have any of the following format (sampleID\_S3\_L001\_R1_001.fastq.gz, sampleID\_1.fastq.gz, sampleID\_1.fq.gz, sampleID\_1.fq .gz), could also be uncompressed.   
* The script `snakenullarbor.sh` creates a new folder where links for the samples are created "option `-o`". In this folder the name of samples is corrected to be read with snakemake in this format "sampleID\_1.fastq.gz"   
* Additionally, `snakenullarbor.sh` writes a config file updated with paths to raw data directory "option `-o`" and snakemake folder "default is the current folder".   
* The option `-fa` refers to additional fasta files for which the fastq files are not avilable but are needed to be included in the analysis. In standard cases, these files should not be compressesd and could end with any of these endings {fasta,fna,fa,fsa,fs,fnn}

#### Example  

`bash ./snakenullarbor.sh -fq data_example/fastq/ -r data_example/Amesancestor.gbk`  
#### Example output of snakenullarbor script 

```
Guessing IDs.....
--------------------------------------------------------------------------------
The follwoing IDs are predicted for the Samples: 
SRR1999417
SRR2968133
SRR2968135


--------------------------------------------------------------------------------
total number of samples to be assembled =3
total number of assemblies =0
--------------------------------------------------------------------------------
Creating links for the fastq reads
ln -s -f /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/data_example/fastq/SRR1999417_1.fastq.gz /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input/SRR1999417_R1.fastq.gz
ln -s -f /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/data_example/fastq/SRR1999417_2.fastq.gz /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input/SRR1999417_R2.fastq.gz
ln -s -f /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/data_example/fastq/SRR2968133_1.fastq.gz /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input/SRR2968133_R1.fastq.gz
ln -s -f /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/data_example/fastq/SRR2968133_2.fastq.gz /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input/SRR2968133_R2.fastq.gz
ln -s -f /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/data_example/fastq/SRR2968135_1.fastq.gz /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input/SRR2968135_R1.fastq.gz
ln -s -f /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/data_example/fastq/SRR2968135_2.fastq.gz /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input/SRR2968135_R2.fastq.gz
--------------------------------------------------------------------------------
Writing reference sequences
seqret -auto -filter -osformat2 fasta < /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/data_example/Amesancestor.gbk > /home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input/ref_strain.fasta
--------------------------------------------------------------------------------
/home/mostafa.abdel/aProjects/gitProjects/snakenullarbor/input is created successfully
--------------------------------------------------------------------------------
writing config file
Found config.yaml. Assume that the current directory is where the snakenullarbor is downloaded. Update the config.yaml file with the paths...
--------------------------------------------------------------------------------
Please note:
To see what snakemake will do, run: snakemake --snakefile snakeNullarbor.Snakefile --cores 128 --use-conda -np 
To execute the pipeline, run: snakemake --snakefile snakeNullarbor.Snakefile --cores 128 --use-conda -p 
To avoid conda problems, run: export PERL5LIB=""
``` 

### 2. Revise and correct the cofig file  
* After running the bash script `snakenullarbor.sh`, the config file "config.yaml" will be updated with paths to raw reads and snakemake_folder
* This file "config.yaml" contains also the settings for different software used in the pipleine. These settings could be adjusted based on the user's preferences 
* __Important__: correct the path to kraken database in the cofig files     



#### 3. Execute the snakenullarbor snakefile 
* The standard output of the `snakenullarbor.sh` involves a hint on how to execute the pipeline with snakemake 

#### Example  
`snakemake --snakefile snakeNullarbor.Snakefile --cores 128 --use-conda -p`

#### Example output report  
Find here [an example report](example/), download the folder locally and open the file index.html using web browser 

## Authors    
Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)  
Jörg Linde (joerg.linde@fli.de)  

## Contributors   
