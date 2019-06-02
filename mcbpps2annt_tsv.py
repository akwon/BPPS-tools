import re
import os
import argparse

#to do -- automate detection of structural location for patterns

parser = argparse.ArgumentParser()
parser.add_argument('-lpr', type = argparse.FileType('r'), required=True) #path and filename for .lpr mcbpps output file
parser.add_argument('-hpt', type=str, required=True) #path and filename for .hpt file used for mcbpps input
parser.add_argument('-struct', type=str, required=True) 
parser.add_argument('-out', type = argparse.FileType('w'), required=True) #path and filename for output tsv file
args = parser.parse_args()

args.out.write('MOTIF_NAME\tCATEGORY#\tSUBGROUP\tOUTGROUP\tPATTERN\tCONSENSUS_POS\tRANK\tFGMATCH\tFGDIVERGE\t%FGMATCH\tBGMATCH\tBGDIVERGE\t%BGMATCH\tINFORMATION_SCORE\tSTRUCTURAL_LOCATION\tKINASE_RESIDUE_FXN\tDESCRIPTION\tEVIDENCE\tREFERENCE\tMOLECULAR_FXN\tREFERENCE(MOL_FXN)\tBIOLOGICAL_FXN\tREFERENCE(BIOL_FXN)\tCOMMENTS\n')

def get_category_name(cat_num):
	with open(args.hpt) as hpt:
		for hptline in hpt.readlines():
			category_line = re.search(cat_num+r'\..*[?!.]',hptline)
			if category_line is not None:
				cat_name = re.search(cat_num+r'\.(.*?)[?!.]',category_line.group(0)).group(1)
				return(cat_name)
				hpt.close()

def get_struct_location(cons_pos):
	with open(args.struct) as structurefile:
		for structline in structurefile.readlines():
			cons_line = re.search(cons_pos+r'\s+.*$',structline)
			if cons_line is not None:
				struct_name = re.search(cons_pos+r'\s+(.*?)$',cons_line.group(0)).group(1)
				return(struct_name)
				structurefile.close()

			
for lprline in args.lpr:   #loops over every line in .lpr file
	category_lpr = re.search(r'BPPS category [0-9]+', lprline) #finds beginning of .lpr output for each category (i.e. :::::::::: BPPS category X: ::::::::::)
	if category_lpr is not None:
		category_num = re.findall(r'[0-9]+',category_lpr.group(0))[0]
		rank = 0
		subgroup_name = get_category_name(category_num)
	else:
		motif_name=outgroup_name=pttrn_aa=pttrn_pos=fg_match=fg_diverge=fgperc=fg_match=bg_diverge=bgperc=info=struct_loc = 'None'
		pttrn = re.search(r'[A-Z]+[0-9]+\([0-9]+\):',lprline) #finds pattern entries for the given category
		if pttrn is not None:
			rank = rank + 1
			#parses columns in .lpr file to find values of variables
			pttrn_aa = re.findall(r'[A-Z]+',pttrn.group(0))[0]
			pttrn_pos = str(re.findall(r'[0-9]+',pttrn.group(0))[0])
			motif_name = str(subgroup_name+'family'+pttrn_aa+pttrn_pos)
			fg_match = str(re.findall(r'[0-9]+',(re.search(r':\s+[0-9]+\s+',lprline)).group(0))[0])
			fg_diverge = str(((re.findall(r'[0-9]+',(re.search(r':\s+[0-9]+\s+[0-9]+\s+\(',lprline).group(0)))))[1])
			fgperc = str((re.findall(r'[0-9]+',((re.findall(r'\([ 0-9]+\%\)',lprline))[0])))[0])
			bg_match =  str(((re.findall(r'[0-9]+',(re.search(r'\([ 0-9]+\%\)\s+[0-9]+\s+',lprline).group(0)))))[1])
			bg_diverge = str(((re.findall(r'[0-9]+',(re.search(r'\([ 0-9]+\%\)\s+[0-9]+\s+[0-9]+',lprline).group(0)))))[2])
			bgperc = str((re.findall(r'[0-9]+',((re.findall(r'\([ 0-9]+\%\)',lprline))[1])))[0])
			info = float(re.findall(r'-?\d+\.\d+',(re.search(r'\([ 0-9]+\%\)\s+[0-9]+\s+[0-9]+\s+\([ 0-9]+\%\)\s+-?\d+\.\d+\(',lprline)).group(0))[0])
			sec_struct = str(get_struct_location(pttrn_pos))
			
			args.out.write(motif_name+'\t'+str(category_num)+'\t'+subgroup_name+'\t'+outgroup_name+'\t'+pttrn_aa+'\t'+pttrn_pos+'\t'+str(rank)+'\t'+fg_match+'\t'+fg_diverge+'\t'+fgperc+'\t'+bg_match+'\t'+bg_diverge+'\t'+bgperc+'\t'+str(info)+'\t'+sec_struct+'\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\n')

			