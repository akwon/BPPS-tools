import os
import subprocess
import glob
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-H', '--hpt', help="hpt file output from omcBPPS (*_new.hpt)", type = argparse.FileType('r'), required=True)
parser.add_argument('-m', '--mma', help="mma file output from omcBPPS (*_new.mma)" , type = str, required=True)
parser.add_argument('-t', '--tweakcma', help="path to tweakcma binary", type = str, required=True)
parser.add_argument('-r', '--rungaps', help="path to run_gaps binary", type = str, required=True)
parser.add_argument('-p', '--profile', help="path to profile used to analyze omcBPPS sets (provide profile name without extension!)", type = str, required=True)
args = parser.parse_args()

write_cmas_command = 'mkdir sets_rungaps; cp ' + args.mma + ' sets_rungaps;' + args.tweakcma + ' sets_rungaps/' + args.mma.replace('.mma','') + ' -write'
subprocess.call(write_cmas_command, shell=True)
for set_cma in glob.glob('sets_rungaps/' + args.mma.replace('.mma','') + '_Set*.cma'):
	write_seq_command = args.tweakcma + ' ' + set_cma.replace('.cma','') + ' -O > ' + set_cma.replace('cma','seq') 
	subprocess.call(write_seq_command, shell=True)
for set_seq in glob.glob('sets_rungaps/' + args.mma.replace('.mma','') + '_Set*.seq'):
	rungaps_command = args.rungaps + ' ' + args.profile + ' ' + set_seq + ' -O'
	subprocess.call(rungaps_command, shell=True)

outfile = open(args.hpt.name + '.annt', 'w')
for hptline in args.hpt:
	if hptline.startswith('+') or hptline.startswith('-') or hptline.startswith('o') and ('-ooooooo') not in hptline:
		set_name = re.search(r'\..*$', hptline).group(0) ; set_name = re.sub(r'[.!?]','', set_name); set_name = re.sub(r'=.*$','', set_name)
		seqcount_list = []
		for cmahit in glob.glob('sets_rungaps/' + args.mma.replace('.mma','') + '_' + set_name + '.seq_*.cma'):
			with open(cmahit) as cmafile:
				cmaline = cmafile.readline().strip()
				if 'aln' in cmahit:
					total_seq_count = re.search(r'\([0-9]+\){', cmaline).group(0).replace(')','').replace('(','').replace('{','')					
				else:
					fam_name = re.search(r'=.*\(', cmaline).group(0).replace('=','').replace('(','')
					fam_seq_count = re.search(r'\([0-9]+\){', cmaline).group(0).replace(')','').replace('(','').replace('{','')
					fam_seq_count = int(fam_seq_count) - 1
					seqcount_list.append(fam_name + ':' + str(fam_seq_count))
					seqcount_list_sorted = sorted(seqcount_list, key=lambda x: int(x.split(':')[1]), reverse=True)
		outfile.write(hptline.strip() + ' ' + total_seq_count + str(seqcount_list_sorted).replace("'", "") + '\n')