import os
from Bio.PDB import *
import argparse
import re

#TO DO: automatically detect pdbstart and pdbend from cma and pdb files.

parser = argparse.ArgumentParser()
parser.add_argument('-pttrn', type = argparse.FileType('r')) #path and filename for .pttrns omcbpps output file
parser.add_argument('-fgname', type = str) #name of BPPS Set. Will be mapped in black->white
parser.add_argument('-category_number', type = str)  #category number of BPPS Set (can be found in .hpt file)
parser.add_argument('-fgname2', type = str) #name of BPPS Set. Will be mapped in dark green->lime green
parser.add_argument('-category_number2', type = str)
parser.add_argument('-fgname3', type = str) #name of BPPS Set. Will be mapped in dark blue->light blue
parser.add_argument('-category_number3', type = str)
parser.add_argument('-fgname4', type = str) #name of BPPS Set. Will be mapped in purple->light pink
parser.add_argument('-category_number4', type = str)
parser.add_argument('-pdb', type = argparse.FileType('r')) #path and filename for pdb structure. structure should be a single chain
parser.add_argument('-cma', type = argparse.FileType('r')) #path and filename to cma aligned pdb sequence
parser.add_argument('-pdbstart', type = int) #residue number in pdb of first residue in the cma alignment
parser.add_argument('-pdbend', type = int) #residue number in pdb of last residue in the cma alignment
parser.add_argument('-pml', type = argparse.FileType('w')) #path and filename for output pml file
args = parser.parse_args()

##Define consensus_array
##The length of consensus_array equals the number of residues in the aligned pdb 
##sequence, and assigns a consensus position number to each residue. If the residue is 
##a consensus-aligned residue, it gets assigned an integer value. If the residue is an
##insert residue, it gets assigned '0'. Gap positions where there is no residue in the
##pdb sequence do not get assigned any value and is not in the array.
for cmaline in args.cma:
	sequence = re.match(r'{\(\)([A-Za-z-]+)\(\)}\*', cmaline)
	if sequence is not None:
		cmaseq = list(sequence.group(1))
		consensus_array = []
		consensus_count = 0
		for aa in cmaseq:
			if aa.isupper():
				consensus_count = consensus_count + 1
				consensus_array.append(consensus_count)
			elif aa == '-':
				consensus_count = consensus_count + 1
			elif aa.islower():
				consensus_array.append(0)
#print(consensus_array)
#print(len(consensus_array))

##Define aligned_res_num_list
parser = PDBParser()
structure = parser.get_structure('struct_name', args.pdb)
residues = structure.get_residues()
res_num_list = []
for res in residues:
	fullid = res.get_full_id()
	res_num = (fullid[3])[1]
	res_num_list.append(res_num)

startindex = res_num_list.index(args.pdbstart)
endindex = res_num_list.index(args.pdbend)
aligned_res_num_list = res_num_list[startindex:endindex+1]
#print(aligned_res_num_list)
#print(len(aligned_res_num_list))

def consensus2pdb(conspos):
	if int(conspos) in consensus_array:
		cons_index = consensus_array.index(int(conspos))
		pdb_num = aligned_res_num_list[cons_index]
	else:
		pdb_num = "none"
	return(pdb_num)

##output a pml file
args.pml.write('# Load the pdb file\n load '+args.pdb.name+'\n hide all\n show cartoon, '+os.path.basename(args.pdb.name).replace('.pdb','')+'\n color white, '+os.path.basename(args.pdb.name).replace('.pdb','')+'\n set cartoon_transparency=0.7, '+os.path.basename(args.pdb.name).replace('.pdb','')+'\n\n')


if args.category_number and args.fgname is not None:
	for line in args.pttrn:
		pttrn_line = re.match(args.category_number+":",line)
		if pttrn_line is not None:
			pttrn_list = (re.sub(args.category_number+": ",'',line.strip('\n'))).split(",")
			selection_string = 'select sel_'+args.fgname+'_cons, '+os.path.basename(args.pdb.name).replace('.pdb','')+' and resi '
			for pttrn in pttrn_list:
				pttrn_aa = re.findall(r'[A-Z]+',pttrn)[0]
				pttrn_consensus_pos = re.findall(r'\d+',pttrn)[0]
				pttrn_rank = pttrn_list.index(pttrn)+1
				pttrn_pdb_pos = consensus2pdb(pttrn_consensus_pos)
				if pttrn_pdb_pos == "none":
					print('Pattern ' + pttrn_aa + pttrn_consensus_pos + ' for ' + args.fgname + ' was not mapped')
					continue

			#write mapped pttrns to pml
				if pttrn_rank <= 10:
					args.pml.write('create '+args.fgname+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname+'_cons'+str(pttrn_rank)+'\ncolor gray10, '+args.fgname+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 10 and pttrn_rank <= 20:
					args.pml.write('create '+args.fgname+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname+'_cons'+str(pttrn_rank)+'\ncolor gray30, '+args.fgname+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname+'_cons'+str(pttrn_rank)+'\n\n')
	
				if pttrn_rank > 20 and pttrn_rank <= 30:
					args.pml.write('create '+args.fgname+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname+'_cons'+str(pttrn_rank)+'\ncolor gray50, '+args.fgname+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 30 and pttrn_rank <= 40:
					args.pml.write('create '+args.fgname+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname+'_cons'+str(pttrn_rank)+'\ncolor gray70, '+args.fgname+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 40 and pttrn_rank <= 50:
					args.pml.write('create '+args.fgname+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname+'_cons'+str(pttrn_rank)+'\ncolor white, '+args.fgname+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname+'_cons'+str(pttrn_rank)+'\n\n')
				
				selection_string = selection_string+str(pttrn_pdb_pos)+'+'
			args.pml.write(selection_string[:-1]+'\n\n')
			args.pttrn.seek(0)
			break
			
if args.category_number2 and args.fgname2 is not None:
	for line in args.pttrn:
		pttrn_line = re.match(args.category_number2+":",line)
		if pttrn_line is not None:
			pttrn_list = (re.sub(args.category_number2+": ",'',line.strip('\n'))).split(",")
			selection_string = 'select sel_'+args.fgname2+'_cons, '+os.path.basename(args.pdb.name).replace('.pdb','')+' and resi '
			for pttrn in pttrn_list:
				pttrn_aa = re.findall(r'[A-Z]+',pttrn)[0]
				pttrn_consensus_pos = re.findall(r'\d+',pttrn)[0]
				pttrn_rank = pttrn_list.index(pttrn)+1
				pttrn_pdb_pos = consensus2pdb(pttrn_consensus_pos)
				if pttrn_pdb_pos == "none":
					print('Pattern ' + pttrn_aa + pttrn_consensus_pos + ' for ' + args.fgname2 + ' was not mapped')
					continue

			#write mapped pttrns to pml
				if pttrn_rank <= 10:
					args.pml.write('create '+args.fgname2+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname2+'_cons'+str(pttrn_rank)+'\ncolor forest, '+args.fgname2+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname2+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 10 and pttrn_rank <= 20:
					args.pml.write('create '+args.fgname2+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname2+'_cons'+str(pttrn_rank)+'\ncolor tv_green, '+args.fgname2+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname2+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 20 and pttrn_rank <= 30:
					args.pml.write('create '+args.fgname2+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname2+'_cons'+str(pttrn_rank)+'\ncolor chartreuse, '+args.fgname2+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname2+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 30 and pttrn_rank <= 40:
					args.pml.write('create '+args.fgname2+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname2+'_cons'+str(pttrn_rank)+'\ncolor limon, '+args.fgname2+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname2+'_cons'+str(pttrn_rank)+'\n\n')
	
				if pttrn_rank > 40 and pttrn_rank <= 50:
					args.pml.write('create '+args.fgname2+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname2+'_cons'+str(pttrn_rank)+'\ncolor paleyellow, '+args.fgname2+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname2+'_cons'+str(pttrn_rank)+'\n\n')
			
				selection_string = selection_string+str(pttrn_pdb_pos)+'+'
			args.pml.write(selection_string[:-1]+'\n\n')
			args.pttrn.seek(0)
			break
	
if args.category_number3 and args.fgname3 is not None:
	for line in args.pttrn:
		pttrn_line = re.match(args.category_number3+":",line)
		if pttrn_line is not None:
			pttrn_list = (re.sub(args.category_number3+": ",'',line.strip('\n'))).split(",")
			selection_string = 'select sel_'+args.fgname3+'_cons, '+os.path.basename(args.pdb.name).replace('.pdb','')+' and resi '
			for pttrn in pttrn_list:
				pttrn_aa = re.findall(r'[A-Z]+',pttrn)[0]
				pttrn_consensus_pos = re.findall(r'\d+',pttrn)[0]
				pttrn_rank = pttrn_list.index(pttrn)+1
				pttrn_pdb_pos = consensus2pdb(pttrn_consensus_pos)
				if pttrn_pdb_pos == "none":
					print('Pattern ' + pttrn_aa + pttrn_consensus_pos + ' for ' + args.fgname3 + ' was not mapped')
					continue

			#write mapped pttrns to pml
				if pttrn_rank <= 10:
					args.pml.write('create '+args.fgname3+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname3+'_cons'+str(pttrn_rank)+'\ncolor density, '+args.fgname3+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname3+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 10 and pttrn_rank <= 20:
					args.pml.write('create '+args.fgname3+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname3+'_cons'+str(pttrn_rank)+'\ncolor blue, '+args.fgname3+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname3+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 20 and pttrn_rank <= 30:
					args.pml.write('create '+args.fgname3+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname3+'_cons'+str(pttrn_rank)+'\ncolor tv_blue, '+args.fgname3+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname3+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 30 and pttrn_rank <= 40:
					args.pml.write('create '+args.fgname3+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname3+'_cons'+str(pttrn_rank)+'\ncolor slate, '+args.fgname3+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname3+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 40 and pttrn_rank <= 50:
					args.pml.write('create '+args.fgname3+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname3+'_cons'+str(pttrn_rank)+'\ncolor lightblue, '+args.fgname3+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname3+'_cons'+str(pttrn_rank)+'\n\n')
				
				selection_string = selection_string+str(pttrn_pdb_pos)+'+'
			args.pml.write(selection_string[:-1]+'\n\n')
			args.pttrn.seek(0)
			break

if args.category_number4 and args.fgname4 is not None:
	for line in args.pttrn:
		pttrn_line = re.match(args.category_number4+":",line)
		if pttrn_line is not None:
			pttrn_list = (re.sub(args.category_number4+": ",'',line.strip('\n'))).split(",")
			selection_string = 'select sel_'+args.fgname4+'_cons, '+os.path.basename(args.pdb.name).replace('.pdb','')+' and resi '
			for pttrn in pttrn_list:
				pttrn_aa = re.findall(r'[A-Z]+',pttrn)[0]
				pttrn_consensus_pos = re.findall(r'\d+',pttrn)[0]
				pttrn_rank = pttrn_list.index(pttrn)+1
				pttrn_pdb_pos = consensus2pdb(pttrn_consensus_pos)
				if pttrn_pdb_pos == "none":
					print('Pattern ' + pttrn_aa + pttrn_consensus_pos + ' for ' + args.fgname4 + ' was not mapped')
					continue
		
			#write mapped pttrns to pml
				if pttrn_rank <= 10:
					args.pml.write('create '+args.fgname4+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname4+'_cons'+str(pttrn_rank)+'\ncolor purpleblue, '+args.fgname4+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname4+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 10 and pttrn_rank <= 20:
					args.pml.write('create '+args.fgname4+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname4+'_cons'+str(pttrn_rank)+'\ncolor deeppurple, '+args.fgname4+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname4+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 20 and pttrn_rank <= 30:
					args.pml.write('create '+args.fgname4+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname4+'_cons'+str(pttrn_rank)+'\ncolor magenta, '+args.fgname4+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname4+'_cons'+str(pttrn_rank)+'\n\n')

				if pttrn_rank > 30 and pttrn_rank <= 40:
					args.pml.write('create '+args.fgname4+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname4+'_cons'+str(pttrn_rank)+'\ncolor violet, '+args.fgname4+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname4+'_cons'+str(pttrn_rank)+'\n\n')
	
				if pttrn_rank > 40 and pttrn_rank <= 50:
					args.pml.write('create '+args.fgname4+'_cons'+str(pttrn_rank)+', resi '+str(pttrn_pdb_pos)+'\n')
					args.pml.write('show sticks, '+args.fgname4+'_cons'+str(pttrn_rank)+'\ncolor lightpink, '+args.fgname4+'_cons'+str(pttrn_rank)+'\nutil.cnc '+args.fgname4+'_cons'+str(pttrn_rank)+'\n\n')
				
				selection_string = selection_string+str(pttrn_pdb_pos)+'+'
			args.pml.write(selection_string[:-1]+'\n\n')
			args.pttrn.seek(0)
			break
