#!/bin/bash -e
#Friedrich-Loeffler-Institut (https://www.fli.de/), IBIZ
#date: April, 2, 2019
#Author: Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)
# A bash script to generate input folder and run the snakefile
#TODO
#directly edit the config file to use with snakemake
#create error if the nunber of IDs less than half the number of the reads, wrong guessing of ids
#call gzip if pigz is not installed, find a way to keep the original files
#if the name of the fasta file are larger than 37 characters, it produces error with prokka
#check dependencies

#dependancies
#(pigz,seqret)
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do cd "$(dirname "$DIR")"; DIR="$(readlink "$(basename "$DIR")")"; done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null #SOURCE
#make a help MSG and pass arguments
PROGNAME=`basename $0`
function usage { echo "USAGE:
   bash ./snakenullarbor.sh -fq fastq_directory -r REF
REQUIRED:
   -fq, --fastq-directory DIR
   -r, --ref 'Reference_strain'
OPTIONAL:
   -fa, --fasta-directory DIR
   -o, --outdir (default: input/)
   -h, --help" ; }
function error_exit { echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2; exit 1; } #SOURCE
function remove_trailing_slash { string="$1" ; new_string=`echo "$string" | perl -nl -e 's/\/+$//;' -e 'print $_'` ; echo $new_string ; } #SOURCE
while [ "$1" != "" ]; do
    case $1 in
        -fq | --fastq-directory ) shift; fastqdirectory=$1 ;;
        -fa | --fasta-directory ) shift; fastadirectory=$1 ;;
        -r | --ref ) shift; reference=$1 ;;
        -o | --outdir ) shift; outdir=$1 ;;
        -h | --help )  usage; exit ;;
        * ) usage; exit 1
    esac
    shift
done
outdir_default='input'
if [[ -z $fastqdirectory ]] && [[ -z $reference ]]; then error_exit "must specify a fastq directory and a reference, use '-fq' and '-r' - exit" ; fi
if [[ -z $fastqdirectory ]] || [[ -z $reference ]]; then error_exit "must specify a fastq directory and a reference, use '-fq' and '-r' - exit" ; fi
if [[ -z $outdir ]]; then outdir=$outdir_default; fi
#variables
fastqdir=`remove_trailing_slash "$fastqdirectory"`
fastadir=`remove_trailing_slash "$fastadirectory"`
outdir=`remove_trailing_slash "$outdir"`
outdir=$( realpath $outdir)
#make the output directory, add a timpestamp if the directoy already exists
#if [[ ! -e $outdir ]]; then mkdir $outdir; fi
#else outdir="$outdir"_"$(date '+%d%b%Y_%H%M%S')" && mkdir $outdir
#fi

mkdir -p $outdir
#first scan the folder for the reads and create error if the reads don't match the required pattern
CheckSamples=$((ls ${fastqdir}/*.{fastq,fastq.gz,fq,fq.gz,fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq )
if [[ $CheckSamples != 1 ]];
  then echo -e "\e[31mcan not guess the samples. Samples names must end with [fastq,fastq.gz,fq,fq.gz] - exit\e[39m";
  exit 1;
fi ;
#get the ID of the sample
printf "Guessing IDs.....\n"
echo "--------------------------------------------------------------------------------"
#id is the what is mentioned before the first underscore in the name
ID=$(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*.{fastq,fq} 2> /dev/null | xargs -n 1 basename 2> /dev/null; ls ${fastqdir}/*{fastq,fq}.gz 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort)
if [[ -e $fastadir ]]; then
  #id is what is mentioned before the last dot for fasta #fasta must be unzipped
  ID_fasta=$(sed 's/\(.*\)\..*/\1/' <(ls ${fastadir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort)
else
  ID_fasta=""
fi
printf "The follwoing IDs are predicted for the Samples: \n${ID}\n${ID_fasta}\n\n"
echo "--------------------------------------------------------------------------------"
printf "total number of samples to be assembled =" && echo ${ID}| wc | awk '{print $2}'
printf "total number of assemblies =" && echo ${ID_fasta}| wc | awk '{print $2}'
echo "--------------------------------------------------------------------------------"
#Compress the fastqfiles if uncompressed
Checkonlyfastq=$((ls ${fastqdir}/*.{fastq,fq} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq )
Checkonlyfq=$((ls ${fastqdir}/*fq 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq )
if [[ $Checkonlyfastq = 1 ]];
  then echo -e "found uncompressed fastq file. Creating .gz files. Original files will not be affected";
  pigz --keep ${fastqdir}/*fastq
  echo "--------------------------------------------------------------------------------"
  if [[ $Checkonlyfq = 1 ]];
    then echo -e "found uncompressed fq file. Creating .gz files. Original files will not be affected";
    pigz --keep ${fastqdir}/*fq
    echo "--------------------------------------------------------------------------------"
  fi;
fi ;
#get the full path of the reads
echo "Creating links for the fastq reads"
for ID in $(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*.{fastq,fq}.gz 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort);
  do
    FILES1=$(realpath $(ls ${fastqdir}/${ID}*_R1_*.gz  2>/dev/null ) 2>/dev/null)
    FILES2=$(realpath $(ls ${fastqdir}/${ID}*_R2_*.gz  2>/dev/null ) 2>/dev/null)
        NF1=$(echo $FILES1 | awk '{print NF}')
    if [[ $NF1 -lt 1 ]]
    then
        #if the reads dont match the format *_R1_*.gz, then check for *_1.fastq. if both formats are not there, then exit
        FILES1=$(realpath $(ls ${fastqdir}/${ID}*_1.*.gz 2>/dev/null ) 2>/dev/null)
        FILES2=$(realpath $(ls ${fastqdir}/${ID}*_2.*.gz 2>/dev/null ) 2>/dev/null)
        NF1=$(echo $FILES1 | awk '{print NF}')
        if [ $NF1 -lt 1 ]
        then
            echo -e "\e[31mfile pattern must match *ID*_R1_*.fastq or *ID*_1.fastq. Files could also be zipped .gz\e[39m"
            exit 1
        fi
    fi
echo "ln -s -f "$FILES1" "$outdir"/"$ID"_R1.fastq.gz"
ln -s -f "$FILES1" "$outdir"/"$ID"_R1.fastq.gz
echo "ln -s -f "$FILES2" "$outdir"/"$ID"_R2.fastq.gz"
ln -s -f "$FILES2" "$outdir"/"$ID"_R2.fastq.gz
done
echo "--------------------------------------------------------------------------------"
#process the already assembled genomes if present
if [[ $((ls ${fastadir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq ) = 1 ]];
  then #printf " \n"; echo "genomes:" >> ${outputfile}; #create the output file
    echo "Creating links for the already assembled genomes"
    for genome in $(sed 's/\(.*\)\..*/\1/' <(ls ${fastadir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort);
      do
        genome_FILE=$(realpath $(ls ${fastadir}/${genome}.{fasta,fna,fa,fsa,fs,fnn}  2>/dev/null ) 2>/dev/null);
        echo "ln -s -f "$genome_FILE" "$outdir"/"$genome".fasta"
        ln -s -f "$genome_FILE" "$outdir"/"$genome".fasta;
      done;
      echo "--------------------------------------------------------------------------------"
  fi
#process the reference strain
if [[ -e $reference ]]; then
  reference=$(realpath $(ls $reference 2>/dev/null ) 2>/dev/null)
  echo "Writing reference sequences"
  echo "seqret -auto -filter -osformat2 fasta < $reference > "$outdir"/Reference.fasta"
  seqret -auto -filter -osformat2 fasta < $reference > "$outdir"/Reference.fasta
fi
#run snakemake command
echo "--------------------------------------------------------------------------------"
echo "$outdir is created successfully"
echo "--------------------------------------------------------------------------------"
#update the config file
echo "writing config file"
if [[ -e config.yaml ]]; then
  config=config.yaml
  echo -e "Found $config. Assume that the current directory is where the snakenullarbor is downloaded. Update the $config file with the paths..."
  var=`pwd`; sed -e "s|snakemake_folder: |snakemake_folder: $var/ #|g" $config > "$outdir"/"$config"
  sed -i "s|raw_data_dir: |raw_data_dir: $outdir/ #|g" "$outdir"/"$config"
  sed -i "s|fasta_dir: |fasta_dir: $outdir/ #|g" "$outdir"/"$config"
  sed -i "s|reference: |reference: $reference #|g" "$outdir"/"$config"
else
  echo -e "\e[1m\e[38:2:240:143:104mCannot find the file: config.yaml. Please update it manually with these paths \nraw_data_dir:\e[0m\e[39m $outdir/ \e[1m\e[38:2:240:143:104m\nfasta_dir:\e[0m\e[39m $outdir/ \e[1m\e[38:2:240:143:104m\nreference:\e[0m\e[39m $reference"
fi

echo "--------------------------------------------------------------------------------"
echo "Please note:"
#echo "Prokka does not like names >37 characters, use 'cd $outdir/ && for f in *.fasta; do echo \${#f}; done' to verify"
echo -e "To see what snakemake will do, run: \e[38;5;42m\e[1msnakemake --snakefile snakeNullarbor.Snakefile --cores 128 --use-conda --keep-going --rerun-incomplete --configfile "$outdir"/"$config" -np \e[39m\e[0m"
echo -e "To execute the pipeline, run: \e[38;5;42m\e[1msnakemake --snakefile snakeNullarbor.Snakefile --cores 128 --use-conda --keep-going --rerun-incomplete --configfile "$outdir"/"$config" -p \e[39m\e[0m"
echo "To avoid conda problems, run: export PERL5LIB=\""\"
#snakemake --snakefile snakeNullarbor.Snakefile -np --cores 128 -p --use-conda
export PERL5LIB=""
