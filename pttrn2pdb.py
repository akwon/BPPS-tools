#!/usr/bin/env python
import os
import subprocess
import sys
from Bio.PDB import *
from Bio import SeqIO 
from Bio.SeqUtils import seq1
import argparse
import re


description = """

***     Requires the following programs and files:       ***
***         BioPython                                    ***
***         MAPGAPS (v2.1.3)                             ***
***         MAPGAPS sequence profile (Default: ePKf)     ***

pttrn2pdb.py: Map residue patterns identified through mcBPPS or omcBPPS onto PDB structures.
  by Annie Kwon
  
  
USAGE: pttrn2pdb.py -p -c -pdb -chain



EXAMPLES:
# Basic
pttrn2pdb.py -p ePK.pttrns -c 1,4 -pdb 1ATP.pdb -chain A -o 1ATP.pml

# Specify set names, colors of residues to be mapped, PDB chains. 
pttrn2pdb.py -p ePK.pttrns -c 1,4 -n ePK,AGC -s blue,purple -pdb 1ATP.pdb -chain A -o 1ATP.pml
""" 
parser = argparse.ArgumentParser(description=description, usage=argparse.SUPPRESS, formatter_class=argparse.RawTextHelpFormatter)
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-p', '--pttrn', type = argparse.FileType('r'), help="path and filename for '.pttrns' file (output from mcBPPS or omcBPPS). This file lists residues unique to each category, ranked by LPR score (in nats). For example, residues unique to category number 2 would be formatted as '2: D175, K30, ... , G10)'", required=True)
requiredNamed.add_argument('-c', '--category', type = str, help="category number(s) of BPPS set to be mapped. Set names and category numbers can be found in '.hpt' file. Multiple categories may be provided in a comma-separated list.", required=True)
requiredNamed.add_argument('-pdb', type = argparse.FileType('r'), help="path and filename for PDB structure.", required=True)
requiredNamed.add_argument('-chain', type = str, help="specify a chain in the PDB structure.", required=True)
parser.add_argument('-o', '--output', type = argparse.FileType('w'), help="path and filename for .pml output file. Default: <pdb file>.pml")
parser.add_argument('-n', '--category_name', type = str, help="name(s) of BPPS set to be mapped. Names should be alphanumeric with words separated by underscores (i.e. EphA3_Eph_TK). As default, category numbers are used as names. Multiple names may be provided in a comma-separated list (should be listed in the same order as the category numbers provided in the '-c' option.") 
parser.add_argument('-s', '--color', type = str, help="colorscale(s) ['blue','purple','green','gray'] of residue sets to be mapped.  Multiple colors may be provided in a comma-separated list (should be listed in the same order as the category numbers provided in the '-c' option.") 
parser.add_argument('-profile', type = str, default="ePKf", help="path and name of sequence profile (should be the same alignment/profile that was used in the mcBPPS/omcBPPS analysis). Uses the ePKf profile by default.") 
parser.add_argument('-cma', type = argparse.FileType('r'), help="path and name of a cma-aligned pdb sequence. The cma alignment should correspond to the same profile used in the BPPS analysis. Providing this skips run_gaps steps.") 
args = parser.parse_args()

###Load the PDB structure and extract the residue sequence and numbering for the specified chain
print("\n\n\n1) Loading the PDB structure and extracting the amino acid sequence...\n...\n...\n...\n")
#def get_pdb_seq(pdbfile):
struct_path = os.path.dirname(args.pdb.name)+'/'
struct_name = os.path.basename(args.pdb.name).replace('.pdb','')
chain_name = os.path.basename(args.pdb.name).replace('.pdb','') +'_'+args.chain.strip() #names the structure and chain
parser = PDBParser(QUIET=True)
structure = parser.get_structure('input_structure', args.pdb) 
chain = structure[0][args.chain]
residues = chain.get_residues()
res_num_list = []; res_seq_list = []

for res in residues:
	fullid = res.get_full_id() 
	res_num = (fullid[3])[1] #parses the residue number 
	res_num_list.append(res_num) 
	res_type = res.get_resname() #gets the residue type in three-letter aa code
	res_seq = seq1(res_type, undef_code='X') #converts to one-letter aa code
	res_seq_list.append(res_seq)
seqfile = open(struct_path+chain_name+'.seq', 'w') #writes PDB sequence to a new file
seqfile.write('>'+chain_name+'\n'+''.join(res_seq_list))
seqfile.close()
#	return(res_num_list, res_seq_list)
#print(res_num_list);print(len(res_num_list));print(res_seq_list);print(len(res_seq_list))

#if a pre-aligned cma is not specified by the -cma option, runs run_gaps on the PDB sequence using ePKf profile (or a profile specified by the '-profile' option)
print("2) Aligning the PDB sequence to the consensus alignment...\n...\n...\n...\n")
if args.cma is None:
	if args.profile is not None:
		rungaps_command = './mapgaps ' + args.profile + ' ' + struct_path+chain_name+'.seq -I=0:0'
	else:
		rungaps_command = './mapgaps ePKf ' + struct_path+chain_name+'.seq + -I=0:0' #needs the '-I' option to not retain flanking segments
	FNULL = open(os.devnull, 'w')
	subprocess.call(rungaps_command, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

#parses the cma file to get the kinase domain start position and aligned kinase domain sequence
if args.cma is not None:
	cmafile = args.cma
else:
	cmafile = open(struct_path+chain_name+'.seq_A.mma', 'r')
for cmaline in cmafile:
	if cmaline.startswith(">"):
		find_startpos = re.search(r'(>.*){\|([0-9]+)(\()', cmaline)
		if find_startpos is None:
			sys.exit("Error: Incorrect cma format. Cannot detect kinase domain start position")
		else:
			startpos = int(find_startpos.group(2))
		seqline = next(cmafile)
		find_sequence = re.match(r'{\(\)([A-Za-z-]+)\(\)}\*', seqline)
		if find_sequence is None:
			sys.exit("Error: Incorrect cma format. Cannot detect aligned kinase domain sequence")		
		else:
			cmaseq = list(find_sequence.group(1))
 
#generates mapping array to convert cma numbering to pdb numbering 		
cma_array = []; consensus_count = 0; pdb_array = []; pdb_count = -1; pdb_seq_array = []
for aa in cmaseq:
	if aa.isupper(): #aligned positions in the aligned pdb sequence
		consensus_count = consensus_count + 1
		pdb_count = pdb_count + 1
		cma_array.append(consensus_count)
		pdb_array.append(res_num_list[startpos+pdb_count])
		pdb_seq_array.append(res_seq_list[startpos+pdb_count])
	elif aa == '-': #gap positions in the aligned pdb sequence
		consensus_count = consensus_count + 1
		cma_array.append(consensus_count)
		pdb_array.append(0)
		pdb_seq_array.append(0)
	elif aa.islower(): #insertion positions in the aligned pdb sequence
		cma_array.append(0)
		pdb_count = pdb_count + 1
		pdb_array.append(res_num_list[startpos+pdb_count])
		pdb_seq_array.append(res_seq_list[startpos+pdb_count])
#print(cma_array);print(len(cma_array));print(pdb_array);print(len(pdb_array))	

def cma2pdb(pos):
	if int(pos) in cma_array:
		pos_index = cma_array.index(int(pos))
		pdb_num = pdb_array[pos_index]
		pdb_res = pdb_seq_array[pos_index]
	else:
		pdb_num = "none"
		pdb_res = "none"
	return(pdb_num,pdb_res)

def get_pttrns(category):
	args.pttrn.seek(0)
	for line in args.pttrn:
		pttrn_line = re.match(category+":",line)
		if pttrn_line is not None:
			pttrn_list = (re.sub(category+": ",'',line.strip())).split(",")
			return(pttrn_list)
	
##parses the .pttrns file, maps patterns to the pdb file, and outputs a pml file
if args.output is not None:
	#outfile = open(args.output, 'w')
	outfile = args.output
else:
	outfile = open(struct_path+chain_name+'.pml', 'w') 
categories = args.category.split(',')
if args.category_name is not None:
	names = args.category_name.split(',')
else:
	names = []
	for cat in categories:
		names.append('Category'+cat)
if args.color is not None:
	colors = args.color.split(',')
	for color in colors:
		if color.lower() not in ['blue','purple','green','gray']:
			print("Error: Incorrect input for -s/--color option.")	
			sys.exit()		
else:
	colors = ['blue','purple','green','gray']*10

def get_colorscale(name,rank):
	if name.lower() == 'blue':
		color_names = ['density','blue','tv_blue','slate','lightblue']
	if name.lower() == 'purple':
		color_names = ['purpleblue','deeppurple','magenta','violet','lightpink']
	if name.lower() == 'green':
		color_names = ['forest','tv_green','chartreuse','limon','paleyellow']
	if name.lower() == 'gray':
		color_names = ['gray10','gray30','gray50','gray70','white']	
	if rank <= 10:
		color_name = color_names[0]
	elif pttrn_rank > 10 and pttrn_rank <= 20:
		color_name = color_names[1]
	elif pttrn_rank > 20 and pttrn_rank <= 30:
		color_name = color_names[2]
	elif pttrn_rank > 30 and pttrn_rank <= 40:
		color_name = color_names[3]
	elif pttrn_rank > 40 :
		color_name = color_names[4]
	return(color_name)

print("3) Mapping BPPS patterns to the PDB structure...\n...\n...\n...\n")
outfile.write('# Load the pdb file\n load '+args.pdb.name+'\n hide all\n show cartoon, '+struct_name+' and chain '+args.chain+'\n color white, '+struct_name+' and chain '+args.chain+'\n set cartoon_transparency=0.7, '+struct_name+'\n\n')
loop_count = -1
for cat in categories:
	loop_count = loop_count + 1
	pttrns = get_pttrns(cat)
	for pttrn in pttrns:
		pttrn_consensus_pos = re.findall(r'\d+',pttrn)[0]
		pttrn_aa = re.findall(r'[A-Z]+',pttrn)[0]
		pttrn_rank = pttrns.index(pttrn)+1
		pttrn_pdb_pos = cma2pdb(pttrn_consensus_pos)		
		#print(pttrn_consensus_pos);print(pttrn_aa);print(pttrn_pdb_pos[0]);print(pttrn_pdb_pos[1])
		if str(pttrn_pdb_pos[0]) == "none":
			print('\tPattern ' + pttrn_aa + pttrn_consensus_pos + ' for ' + names[loop_count] + ' was not mapped')
			continue
		if str(pttrn_pdb_pos[0]) == "0":
			print('\tPattern ' + pttrn_aa + pttrn_consensus_pos + ' for ' + names[loop_count]  + ' is not ordered in the provided PDB structure')
			continue
		if str(pttrn_pdb_pos[1]) not in list(pttrn_aa):
			print('\tPattern ' + pttrn_aa + pttrn_consensus_pos + ' for ' + names[loop_count]  + " was found as a --'" + str(pttrn_pdb_pos[1]) +"'-- in the provided PDB structure")
		#write mapped pttrns to pml
		mapping_color = get_colorscale(colors[loop_count], pttrn_rank)
		outfile.write('create '+names[loop_count]+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos[0])+'\n')
		outfile.write('show sticks, '+names[loop_count]+'_cons'+str(pttrn_rank)+'\ncolor ' + mapping_color + ', ' + names[loop_count]+'_cons'+str(pttrn_rank)+'\nutil.cnc '+names[loop_count]+'_cons'+str(pttrn_rank)+'\n\n')
outfile.close()

print("\n\n\nFINISHED! Wrote output file to: " + str(os.path.join(outfile.name)) + "\n\n\n\n")
