#!/usr/bin/env python
import os
import subprocess
import glob
import argparse
import re
import numpy as np
import sys
from Bio import Entrez
from ete3 import NCBITaxa

description = """

***     Requires the following programs and files:       ***
***         ete3/NCBI taxonomy database                  ***
***         tweakcma                                     ***
***         MAPGAPS (v2.1.3)                             ***
***         MAPGAPS sequence profile (Default: ePKf)     ***

annt_parti.py: Annotate sequence partitions identified through mcBPPS or omcBPPS. The input for this program is a partitioned alignment file that is generated as mcBPPS/omcBPPS output ('_new.mma'). There should be a corresponding '_new.hpt' file in the same directory.
  by Annie Kwon
  
  
USAGE: annt_hpt.py <input>

     <input>     path and file name to a '_new.mma' file


EXAMPLES:
# Basic
annt_parti.py <input>

# Specify a different profile
annt_parti.py <input> -p TK-Group

# Print a list of species represented in each partition and print a species tree
annt_parti.py <input> -s -tree

IF RUNNING ON SAPELO:
    module load Biopython
    module load ete3
    module load numpy
""" 
parser = argparse.ArgumentParser(description=description, usage=argparse.SUPPRESS, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', '--output', help="output annotated .hpt file", type = str)
parser.add_argument('-tsv', help="output a .tsv file with annotated partition information", type = str)
parser.add_argument('-p', '--profile', help="path to MAPGAPS sequence profile (provide profile name without extension)", type = str, default = './ePKf')
parser.add_argument('-t', '--tweakcma', help="path to tweakcma binary", type = str, default = './tweakcma')
parser.add_argument('-m', '--mapgaps', help="path to mapgaps binary", type = str, default = './mapgaps')
parser.add_argument('-s', '--species', help="print the list of species represented in each partition",action="store_true")
parser.add_argument('-tree', help="print the species tree of species represented in each partition",action="store_true")

args, unk = parser.parse_known_args()
if len(unk) != 1:
    sys.exit('Input file not provided')
input_mma = unk[0]
mma_path = os.path.dirname(input_mma)+'/'
mma_name = os.path.basename(input_mma).replace('.mma','')
new_directory = mma_path + 'partition_cmas/'
if args.output is not None:
	outfile = open(args.output, 'w')
else:
	outfile = open(mma_path + mma_name + '.hpt.annt' , 'w')
if args.tsv is not None:
	tsvfile = open(args.tsv, 'w') 


def write_cmas(mma_file): #divides partitions in '_new.mma ' file into individual '.cma' files
	write_cmas_command = 'mkdir ' + new_directory + '; cp ' + mma_file + ' ' + new_directory + ';' + args.tweakcma + ' ' + new_directory + mma_name + ' -write'
	FNULL = open(os.devnull, 'w')
	subprocess.call(write_cmas_command, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

def write_seqs(cma_file):	
	write_seq_command = args.tweakcma + ' ' + cma_file.replace('.cma','') + ' -O > ' + cma_file.replace('.cma','.seq') 
	FNULL = open(os.devnull, 'w')
	subprocess.call(write_seq_command, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

def get_cmainfo(cma_file):
	count_regex = r'\([0-9]+\){'
	name_regex = r'=.*\('
	with open(cma_file) as cma:
		cmaline = cma.readline().strip()
		fam_name = re.search(name_regex, cmaline).group(0).replace('=','').replace('(','')
		seq_count = re.search(count_regex, cmaline).group(0).replace(')','').replace('(','').replace('{','')
	return(seq_count,fam_name)

def get_tax(cma_file): #finds the 'lowest common ancestor' of species represented in a cma file  
	ncbi = NCBITaxa()
	org_regex = r'\[(.*)\]'
	taxid_list = []
	for line in open(cma_file,'r'):
		if line.startswith(">"):
			find_org_name = re.search(org_regex,line)
			if find_org_name is not None:
				org_name = find_org_name.group(1)
				taxid = str(ncbi.get_name_translator([org_name])); taxid = re.sub(r'^.*\[','',taxid); taxid = re.sub(r'\].*$','',taxid)
				if taxid != '{}' and taxid != '32630' and taxid != '10239': #omit sequences from viruses and synthetic constructs'
					taxid_list.append(taxid)
	tax_list = ncbi.get_taxid_translator(taxid_list)
	tree = ncbi.get_topology(taxid_list) 
	tree_labeled = tree.get_ascii(attributes=['sci_name', 'taxid'])
	lca_id = str(tree.get_tree_root) ; lca_id = re.sub(r"^.*node '",'',lca_id) ; lca_id = re.sub(r"'.*$",'',lca_id)
	lca_name = str(ncbi.get_taxid_translator([lca_id])); lca_name = re.sub(r"'}$",'',lca_name) ; lca_name = re.sub(r"^.*'",'',lca_name)
	return(lca_name,tax_list,tree_labeled)
	
def align_seqs(seq_file):
	mapgaps_command = args.mapgaps + ' ' + args.profile + ' ' + seq_file + ' -O'
	FNULL = open(os.devnull, 'w')
	subprocess.call(mapgaps_command, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

def get_cmahits(set): #makes a sorted list of profile hits for a given partition (sorted by the # of sequence hits)
	total_seq_count = get_cmainfo(new_directory + mma_name + '_' + set + '.seq_A.mma')[0]
	seqcount_list = []
	for cmahit in glob.glob(new_directory + mma_name + '_' + set + '.seq_*.cma'):
		fam_seq_count = int(get_cmainfo(cmahit)[0]) - 1 #subtract 1 to account for the consensus sequence that gets added via MAPGAPS
		fam_name = get_cmainfo(cmahit)[1]
		seqcount_list.append(fam_name + ':' + str(fam_seq_count))
	seqcount_list_sorted = sorted(seqcount_list, key=lambda x: int(x.split(':')[1]), reverse=True)
	return(total_seq_count, seqcount_list_sorted)

print("\n\n\n1) Extracting sequences from each partition...\n...\n...\n...\n")	
write_cmas(input_mma)
for set_cma in glob.glob(new_directory + mma_name + '_*.cma'):
	write_seqs(set_cma)	

print("2) Comparing each partition to families in the sequence profile " + args.profile + "...\n...\n...\n...\n")
for set_seq in glob.glob(new_directory + mma_name + '_*.seq'):
	align_seqs(set_seq)

print("3) Gathering taxonomic information for each partition...\n...\n...\n...\n")	
hptfile = open(mma_path + mma_name + '.hpt' , 'r')
for hptline in hptfile:
	if hptline.startswith('+') or hptline.startswith('-') or hptline.startswith('o'):
		set_name = re.search(r'\..*$', hptline).group(0) ; set_name = re.sub(r'[.!?]','', set_name); set_name = re.sub(r'=.*$','', set_name)
		if set_name != 'Random':
			annt = get_cmahits(set_name)
			tax = get_tax(new_directory + mma_name + '_' + set_name + '.cma')
			outfile.write(hptline.strip() + ' ' + annt[0] + str(annt[1]).replace("'", "") + '\t' + str(tax[0]) + '\n')
			if args.species == True:
				species_list = []
				#for species in list(str(tax[1])):
				for species in str(tax[1]).split(','):
					species = re.sub(r"^.*u'",'',species); species = re.sub(r"'.*$",'',species) #do some processing to only get the species name
					species_list.append(species)
				print (set_name + '\t' + str(species_list) + '\n')
			if args.tree == True:
				print (set_name + '\n' + str(tax[2]))
			if args.tsv is not None:
				tsvfile.write(set_name + '\t' + annt[0] + '\t' + str(annt[1]).replace("'", "").replace("[","").replace("]","") + '\t' + str(tax[0]) + '\n')
		
print("4) Wrote output file: " + outfile.name + ".\n\n\n")