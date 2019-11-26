#!/bin/bash -e
#tmp_file=$mktemp
file=$1
#out=$2
cat $file  | sed -e "s/[][,']//g;s/[)()]//g" | tr ' ' '\n' | sed '/ref_strain/d' #> $out

#cat $file  | sed -e "s/[][,']//g;s/[)()]//g" | tr ' ' '\n' #> $out
