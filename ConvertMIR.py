#!/usr/bin/env python

# Sarah B. Kingan
# 27 May 2015
# 
# miRNA project
# 
# This python script makes species-specific miRNA sequence libraries
# Input: miRNA sequences from D. melanogaster in miREexpress format
# Output: miRNA sequences in miRExpress format for D. simulans clade species

###########################################################################################################
# import libraries
from argparse import ArgumentParser
from subprocess import check_output, call
import csv
import string
import sys

###########################################################################################################
# define command line arguments, program description, and help
desc ='''Convert D.mel miRNA sequences into species-specific sequence
		by aligning D.mel miRNA sequences to miRNA precursor sequences from simulans clade.'''
parser = ArgumentParser(description=desc)
parser.add_argument("miRNA_precursorFilePath", help="path to species-specific precursor sequences")
parser.add_argument("speciesID", help="species identifier, e.g. dmau, dsim, dsec")
args = parser.parse_args()

###########################################################################################################
# get filenames for blast database, genome sequence in fasta format and D.melanogaster miRNA precursor seqs
miRNA_precursorFile = args.miRNA_precursorFilePath
species = args.speciesID
dmelMIRFile = "/home/LCPG/skingan/small_RNA/miRExpress/fourSpLibraries/final_files/dmel_miRNA.txt"
open(species + '_miRNA.txt', 'w')
open(species + '_miRNA.err', 'w')




###########################################################################################################
# Program Body	
###########################################################################################################

### blast (bl2seq) each Dmel miRNA seq to corresponding sim-clade precursor sequence ###

# search for precursor seq in species file	
with open(dmelMIRFile, "r") as mirFile:
	for mirLine in csv.reader(mirFile, delimiter='\t'):
		seq = ''
		miR_name = mirLine[0][:-3]
		with open(miRNA_precursorFile, "r") as precFile:
			for precLine in precFile:
				miR_search = miR_name.replace('R','r') + "\t"
				if miR_search in precLine:
					precLineList = precLine.split("\t")
# write tmp miRNA fasta file for blast
					with open('miRNA.fa', 'w') as mirFasta:
						mirFasta.write('>' + mirLine[0] + "\n") 	# miRNA name
						mirFasta.write(mirLine[1] + "\n") 			# miRNA seq
# write tmp precursor file for blast
					with open('precursor.fa', 'w') as precFasta:
						precFasta.write('>' + precLineList[0] + "\n") 	# precursor name
						precFasta.write(precLineList[1]) 				# precursor seq
# perform blast between mel miRNA and sim-clade precursor IF precursor seq is present
					command = 'bl2seq -i miRNA.fa -j precursor.fa -p blastn -D 1'
					bo = check_output(command, shell=True)
# if blast returns something, process best hit
					boList = bo.split("\n")
					if len(boList) > 4: # there is a blast hit, 3 header lines, 1 results line, 1 blank line
						hit = boList[3] # best blast hit
						hitList = hit.split("\t")
# if best blast hit is full length and 100% ID, print mel miRNA
						miR_length = len(mirLine[1])
						if int(hitList[8]) < int(hitList[9]):
							if int(miR_length) == int(hitList[3]) and hitList[2] == '100.00':
								seq = mirLine[1]
								start_pos = hitList[8]
# if blast hit has mismatch, slice sim-clade precursor seq and print to file
							else:
								if float(hitList[3])/float(miR_length) > 0.8 and int(hitList[8]) < int(hitList[9]):
									start = max(int(hitList[8]) - int(hitList[6]),0)
									end = int(hitList[9]) + (int(hitList[7]) - int(miR_length))
									seq = precLineList[1][start:end]
									start_pos = start + 1
							if len(seq) > 0:
								with open(species + '_miRNA.txt', 'a') as out:
									out.write(mirLine[0] + "\t" + seq + "\t" + precLineList[0] + "\t" + str(start_pos) + "\n")
						else:
							with open(species + '_miRNA.err', 'a') as err:
								err.write(mirLine[0] + "\t" + precLineList[0] + "\n")
					else:
						with open(species + '_miRNA.err', 'a') as err:
							err.write(mirLine[0] + "\t" + precLineList[0] + "\n")
							
							


