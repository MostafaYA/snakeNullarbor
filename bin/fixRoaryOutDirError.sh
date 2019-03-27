#!/bin/bash
#Mostafa.Abdel-Glil@fli.de, date: October, 24, 2018
# check if a directory exists and move it to another directory
# the script was made for a specific use within the pangenome.Snakemake workflow. Thus no time stamp subfolder will be created for Roary

directory=$1
outdirectory=$2
if [ -d $directory ] ; then echo "directory exists" ; mkdir -p $outdirectory; mv  $directory $outdirectory; fi
