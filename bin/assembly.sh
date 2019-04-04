#!/bin/bash -e
#Friedrich-Loeffler-Institut (https://www.fli.de/)
#date: 21.03.2019
#Author: Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)
#dependencies -> #shovill, spades, trimmomatic
#TODO (MaSuRCA, A5, velvet, correct the assembly using pilon, reference based assembly, quast)
#get to the script directory
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do cd "$(dirname "$DIR")"; DIR="$(readlink "$(basename "$DIR")")" ; done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null #SOURCE

PROGNAME=`basename $0`
function usage { echo "USAGE:
  bash ./assembly.sh --R1 reads1.fastq.gz --R2 reads2.fastq.gz
REQUIRED:
 --R1 read1
 --R2 read2
OPTIONAL:
  -a, --assembler [spades shovill megahit sksea] (default: shovill)
  -d, --outdir (default: '\$assembler')
  -o, --output (default: contigs.fa)
  -p, --threads (default: 16)
  -t, --trim [trimmomatic]
  -f, --filter [true] (length > 500 & kmer coverage > 3)
  -opts, --assembler_option
  -h, --help " ; }
function error_exit { echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2; exit 1; } #SOURCE
function remove_trailing_slash { string="$1" ; new_string=`echo "$string" | perl -nl -e 's/\/+$//;' -e 'print $_'` ; echo $new_string ; } #SOURCE

#defaults
assembler_default="shovill"
#OUTDIR_default=$assembler_default , assembler name
CONTIGSoutput_default='contigs.fa'
THREADS_default='16'
#variables
while [ "$1" != "" ]; do
    case $1 in
        -a | --assembler ) shift; assembler=$1 ;;
        --R1 ) shift; READ1=$1 ;;
        --R2 ) shift; READ2=$1 ;;
        -d | --outdir ) shift; OUTDIR=$1 ;;
        -o | --output ) shift; CONTIGSoutput=$1;;
        -p | --threads ) shift; THREADS=$1 ;;
        -t | --trim ) shift; TRIM=$1 ;;
        -f | --filter ) shift; FILTER=$1 ;;
        -opts | --assembler_option )  shift; assembler_options=${1} ;;
        -h | --help )  usage; exit ;;
        * ) usage; exit 1
    esac
    shift
done
#if [ $# -lt 1 ]; then usage; exit  ; fi
if [[ -z $READ1 ]] && [[ -z $READ2 ]] ; then error_exit "Reads are required. use --R1 and --R2" ; fi
if [[ -z $READ1 ]] | [[ -z $READ2 ]] ; then error_exit "Reads are required. use --R1 and --R2" ; fi
if [[ -z $assembler ]]; then assembler=$assembler_default; fi
if [[ -z $OUTDIR ]]; then OUTDIR=$assembler; fi
if [[ -z $CONTIGSoutput ]]; then CONTIGSoutput=$CONTIGSoutput_default; fi
if [[ -z $THREADS ]]; then THREADS=$THREADS_default; fi
if ! [[  $THREADS =~ '^[0-9]+$' ]]; then THREADS=$THREADS_default; fi

#processing variables
PROCDIR=$(pwd)
temp_dir=$(mktemp -d)
OUTDIR=`remove_trailing_slash "$OUTDIR"`
DIR=`remove_trailing_slash "$DIR"`
LIBDIR="$DIR/../lib"

#assembler=shovill
if [[ "$assembler" == "shovill" ]];
then
    echo "using shovill for assembly"
    if [[ "$TRIM" = "trimmomatic" ]] && [[  $FILTER =  "true" ]];
    then
      assembler_options_additional='--trim --minlen 500 --mincov 3'
    elif [[ "$TRIM" = "trimmomatic" ]] && [[  "$FILTER" !=  "true" ]];
    then
      assembler_options_additional='--trim'
    elif [[ "$TRIM" != "trimmomatic" ]] && [[  "$FILTER" =  "true" ]];
    then
      assembler_options_additional='--minlen 500 --mincov 3'
    else
      assembler_options_additional=""
    fi
    echo "shell command:  shovill --outdir $OUTDIR --cpus $THREADS --R1 $READ1 --R2 $READ2 --force $assembler_options $assembler_options_additional"
    shovill --outdir $OUTDIR --cpus $THREADS --R1 $READ1 --R2 $READ2 --force $assembler_options $assembler_options_additional
    cp -r $OUTDIR/contigs.fa ./$CONTIGSoutput
fi

#assembler=spades
#spades/assembly_graph_with_scaffolds.gfa' => 'contigs.gfa'
if [[ "$assembler" == "spades" ]];
then
    echo "using SPAdes for assembly"
    #OUTDIR=$assembler #
    if [[ "$TRIM" = "trimmomatic" ]];
    then
      echo "trimming using trimmomatic"
      echo "shell command: trimmomatic PE -threads $THREADS -phred33 $READ1 $READ1  $OUTDIR/R1.fq.gz $OUTDIR/R2.fq.gz  LEADING:10 TRAILING:10 MINLEN:30 ILLUMINACLIP:$LIBDIR/trimmomatic.fa:1:30:11"
      mkdir -p $OUTDIR &&
      trimmomatic PE -threads $THREADS -phred33 $READ1 $READ1  $OUTDIR/R1.fq.gz /dev/null $OUTDIR/R2.fq.gz  /dev/null ILLUMINACLIP:$LIBDIR/trimmomatic.fa:1:30:11 LEADING:10 TRAILING:10 MINLEN:30 2>&1 | sed 's/^/[trimmomatic] /' &&
      READ1=$OUTDIR/R1.fq.gz
      READ2=$OUTDIR/R2.fq.gz
    fi
    echo "shell command: spades.py --pe1-1 $READ1 --pe1-2 $READ2 --threads $THREADS  --careful --memory 32 -o $temp_dir --tmp-dir $(mktemp -d) -k 31,55,79,103,127 $assembler_options "
    spades.py	--pe1-1	$READ1	--pe1-2 $READ2	--threads	$THREADS --careful --memory	32	-o	$temp_dir	--tmp-dir	$(mktemp -d)	-k	31,55,79,103,127 $assembler_options 2>&1 | sed 's/^/[spades] /' &&
    mkdir -p $OUTDIR
    if [[  $FILTER =  "true" ]];
    then
      python3.5 $DIR/filter_contigs.py -i $temp_dir/contigs.fasta -o $CONTIGSoutput -c 3 -l 500
    else
        cp -v -f "$temp_dir/contigs.fasta" "$CONTIGSoutput"
    fi
    cp -v -f "$temp_dir/"{spades.log,scaffolds.fasta,contigs.fasta} "$OUTDIR/"
    cp -v -f "$temp_dir/assembly_graph_with_scaffolds.gfa" "$OUTDIR/scaffolds.gfa"
    cp -v -f "$temp_dir/assembly_graph.fastg" "$OUTDIR/contigs.gfa"
    rm -frv "$temp_dir"
fi



#READ1=snakemake@input[["r1"]]
#READ2=snakemake@input[["r2"]]
#OUTDIR=snakemake@params[["assembly_dir"]]
#CONTIGSoutput=snakemake@output[["contig"]]
#THREADS=snakemake@threads[[1]]
