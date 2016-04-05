#!/usr/bin/env python

# Sarah B. Kingan
# 2 June 2015
# 
# miRNA project
# 
# This python script corrects the starting files for dme MIR sequence from miRExpress
# Input: dme_precursor.txt dme_miRNA.txt
# Output: corrected versions of files

###########################################################################################################
# import libraries
from argparse import ArgumentParser
from subprocess import check_output, call
import csv
import string
import sys


###########################################################################################################
# get filenames
precursorFilename = 'dme_precursor.txt'
miRNAFilename = "dme_miRNA.txt"


###########################################################################################################
# Program Body	
###########################################################################################################


# search for precursor seq in species file	
with open(miRNAFilename, "r") as miRNAFile:
	keys = ['miRNA', 'seq', 'gene', 'position']
	dict = {}
	new_key = ''
	for line in csv.reader(miRNAFile, delimiter='\t'):
		child = {}
		for k,v in zip(keys,line):
			child[k] = v
		new_key = line[2].replace('r','R') + line[0][-3:]	
		dict[new_key] = child
#print dict		

for miRNA in sorted(dict.keys()):
	print miRNA + "\t" + dict[miRNA]['seq'] + "\t" + dict[miRNA]['gene'] + "\t" + dict[miRNA]['position']

	


