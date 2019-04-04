#!/bin/bash
#Mostafa.Abdel-Glil@fli.de, date 2018.10.22
#script to replace the fasta header with the file name prefix and sequential numbers

file=$1
output=$2

#bioawk -v filename=$(awk 'BEGIN{FS="."}{ print $1 }' <(basename $file)) -c fastx '{ print ">"filename"\n"$seq; }' $file |\
# awk '/^>/{$0=$0"."(++i)}1' > $output
 bioawk -v filename=contig -c fastx '{ print ">"filename"\n"$seq; }' $file |\
  awk '/^>/{$0=$0"."(++i)}1' > $output
