#!/bin/bash -e
#Friedrich-Loeffler-Institut (https://www.fli.de/)
#date: 21.03.2019
#Author: Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)
#dependencies -> #shovill, spades, trimmomatic

#get to the script directory
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

PROGNAME=`basename $0`
function usage { echo "USAGE:
  bash ./assembly.sh --R1 reads1.fastq.gz --R2 reads2.fastq.gz
REQUIRED:
 --R1 read1
 --R2 read2
 --db database
OPTIONAL:
  -t, --taxoner [kraken kraken2 centrifuge metaphlan]
  -o, --output (kraken.tab)
  -p, --threads
  -opts, --taxoner_option
  -h, --help " ; }
function error_exit { echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2; exit 1; }
function remove_trailing_slash { string="$1" ; new_string=`echo "$string" | perl -nl -e 's/\/+$//;' -e 'print $_'` ; echo $new_string ; }

#defaults
taxoner_default="kraken"
taxonomyOutput_default='kraken.tab'
THREADS_default='16'


#processing variables
PROCDIR=$(pwd)
DATABASE=
OUTDIR=`remove_trailing_slash "$OUTDIR"`
DIR=`remove_trailing_slash "$DIR"`
LIBDIR="$DIR/../lib"

#taxoner=kraken
if [[ "$taxoner" == "kraken" ]];
then
    echo "shell command:  kraken --db $DATABASE --paired --check-names --threads $THREADS --gzip-compressed --fastq-input $READ1 $READ2 | kraken-report > $taxonomyOutput"
    kraken --paired --db $DATABASE --check-names --threads $THREADS --gzip-compressed --fastq-input $READ1 $READ2 | kraken-report >  "$taxonomyOutput"
fi
