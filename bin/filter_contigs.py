#!/usr/bin/python
# filter SPAdes contigs based on header for coverage and minimum length  
# to obtain nucleotide sequences for all genes 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os.path

def handle_args():
    import argparse
    usage = ""
    usage += "Filter SPAdes contigs for coverage and minimum length"

    parser = argparse.ArgumentParser(description = usage )
    #parser.print_help()
    parser.add_argument('-i', dest='infile', help = 'fasta [*.fasta]')
    parser.add_argument('-o', dest='outfile', help='FASTA file [*.fasta]')
    parser.add_argument('-c', dest='min_coverage', type=float, help='minimum coverage [float]')
    parser.add_argument('-l', dest='min_length', type=int, help='minimum length [int]')
    args = parser.parse_args(sys.argv[1:])
  #  args = parser.parse_args()
    return args



def main(options):
    # parse multiple geneBank records (as a single SeqRecord object),
    # each with multiple features (as SeqFeature objects)
    input_filename = options.infile
    output_filename = options.outfile
    input_handle = open(input_filename, "r")
    output_handle = open(output_filename, "w")
    
    
    filename = os.path.basename(output_filename)
    records = list(SeqIO.parse(input_handle, "fasta"))
    seqlen = 0
    num_filtered = 0
    
    for i in range(len(records)):
      seqlen += len(records[i].seq)
      IDa = records[i].id.split('_')
      contigID = IDa[0] + "_" + IDa[1]
      contigLength = int(IDa[3])
      contigCoverage = float(IDa[5])
      if len(records[i].seq) >= options.min_length and contigCoverage >= options.min_coverage:
          SeqIO.write(records[i], output_handle, "fasta")
          # print records[i].id
          # print contigID + "\t" + str(len(records[i].seq)) + "\t" + str(contigCoverage)
          # print "len:" + str(len(records[i].seq))
          # print "length:" + str(contigLength)
          # print "coverage:" + str(contigCoverage)
          num_filtered += 1
    
    print ("processed " + str(len(records)) + " records, total length " + str(seqlen) + "\n")
    print ("filtered " + str(num_filtered) + " contigs that are at least " + str(options.min_length) + " bp long and with coverage >= " + str(options.min_coverage))
    
    output_handle.close()     
    input_handle.close()

if __name__ == '__main__':
    options = handle_args()
    main(options)

