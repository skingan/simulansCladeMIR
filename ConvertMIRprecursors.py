#!/usr/bin/env python

# Sarah B. Kingan
# 27 May 2015
# 
# miRNA project
# 
# This python script makes species-specific miRNA precursor sequence libraries
# Input: miRNA precursor sequences from D. melanogaster from miRExpress program
# Output: miRNA precursor sequences in miRExpress format for user-specified species

###########################################################################################################
# import libraries
from argparse import ArgumentParser
from subprocess import check_output
import csv
import string
import sys

###########################################################################################################
# define command line arguments, program description, and help
desc ='''Convert D.mel miRNA sequences into species-specific sequence
		by BLASTing D.mel sequences to genome sequences of the simulans clade.'''
parser = ArgumentParser(description=desc)
parser.add_argument("blastDB", help="path to blastDB (basename)")
parser.add_argument("genomeSeq", help="path to genome sequence in fasta format")
parser.add_argument("speciesID", help="species identifier, e.g. dmau, dsim, dsec")
args = parser.parse_args()

###########################################################################################################
# get filenames for blast database, genome sequence in fasta format and D.melanogaster miRNA precursor seqs
DB = args.blastDB
genome = args.genomeSeq
species = args.speciesID
dmelMIR = "dme_precursor.txt"
open(species + '_precursor.txt', 'w')
open(species + '_precursor.err', 'w')


###########################################################################################################
# input blast output and return coordinates of best hit (or note)
def blast2coords(blastOut):
	if len(blastOut) < 11: # there is no blast hit
		return str('no hit')
	elif float(blastOut[3])/float(blastOut[2]) >= 0.9: # alignment is >= 90% miRNA length
# get chrom and coordinates
		chrom = blastOut[1]
		if float(blastOut[6]) < float(blastOut[7]): 
			start = float(blastOut[6]) - (float(blastOut[4]) - 1)
			end = float(blastOut[7]) + (float(blastOut[5]) - float(blastOut[2]))
		else:
			start = float(blastOut[6]) + (float(blastOut[4]) - 1)
			end = float(blastOut[7]) - (float(blastOut[5]) - float(blastOut[2]))
		return [chrom,int(start),int(end)]
	else: # alignment is too short
		return str('alignment too short')

###########################################################################################################
# convert coords to sequence using samtools faidx
def coords2seq(fastaFile,chrom,start,end):
	command = 'samtools faidx ' + fastaFile + ' ' + chrom + ':' + str(min(start,end)) + '-' + str(max(start,end))
# get output from samtools faidx command
	so = check_output(command, shell=True)
# format sequence header and RNA message sequence
	seqOutput = so.split("\n")
	h = seqOutput[0]
	header = h.replace(">","")
	s = seqOutput[1:]
	seq = ''.join(s)
	RNA = seq.translate(string.maketrans("T","U"))
	if start > end:
		R = ''.join(reversed(RNA))
		RC = R.translate(string.maketrans("ACGU","UGCA"))
		RNA = RC
	return[header,RNA]


###########################################################################################################
# Program Body	
###########################################################################################################

# blast each Dmel miRNA seq to specified database
with open(dmelMIR) as f:
	for l in csv.reader(f, delimiter='\t'):
# write tmp fasta file
		with open('tmp.fa', 'w') as t:
			t.write('>' + l[0])
			t.write("\n")
			t.write(l[1])
			t.write("\n")
# blast tmp file to specified database and save output (str with line break characters)
		bo = check_output('blastn -max_target_seqs 1 -outfmt \"6 qseqid sseqid qlen length qstart qend sstart send pident gaps evalue\" -query tmp.fa -db ' + DB, shell=True) 
		blastOutput = bo.split("\t")
# extract chromosome coordinates
		coords = blast2coords(blastOutput)
# print sequences (UC) to .txt file
		if type(coords) is list: 
			hs = coords2seq(genome,coords[0],coords[1],coords[2])
			U = hs[1].upper()
			with open(species + '_precursor.txt', 'a') as out:
				out.write(l[0] + "\t" + U + "\n")
# print miRNA without blast hits to .err file
		else:
			with open(species + '_precursor.err', 'a') as err:
				err.write(l[0] + "\t" + coords + "\n")


