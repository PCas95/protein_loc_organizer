#!/usr/bin/env python

# tmhmm_extract.py
__version__ = '1.26.01'

# manage libraries
import argparse, sys, os, re
from datetime import datetime
from collections import defaultdict
from common_utils import id_reader, find_id


# classes
# custom error classes
class UnreadableInputError(ValueError):
	"""Input file ureadable or empty"""

class SeqIDNotFoundError(ValueError):
	"""No matching SeqIDs found in input"""


# globals
silent = False
header_map = {'len': 'Length',
			'PredHel': 'Number of predicted TMHs',
			'ExpAA': 'Exp number of AAs in TMHs',
			'First60': 'Exp number, first 60 AAs'}


# functions
## saves output table in standardised format
def save_output(ofi: str, ret_obj: dict[list]):
	## get full column range for output table
	nCol_max = 0
	id_max_col = ''
	for x in ret_obj.keys():
		current_nCol = len(ret_obj[x].keys())
		if current_nCol > nCol_max:
			nCol_max = current_nCol
			id_max_col = x
	## write table
	with open(ofi, 'w') as out:
		print('SeqID', "\t".join(list(ret_obj[id_max_col].keys())), sep="\t", file=out)
		for c in ret_obj.keys():
			values = [ret_obj[c][x] for x in list(ret_obj[id_max_col].keys())]
			print(c, *values, sep="\t", file=out)


## warning/error check function
def checkers(dc: dict):
	# produce warning log + message if there are no results for some SeqIDs
	errList = []
	for k in dc.keys():
		if len(dc[k]) < 6:
			errList.append(str(dc[k]))
	if len(errList) == len(dc.keys()):
		print('ERROR: It seems that none of the provided SeqIDs has TMHMM prediction metrics.\nPlease double check the input file.')
		sys.exit(1)
	elif errList == True:
		warnings = f'warnings_tmhmm_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.txt'
		print(f'WARNING: One or more SeqIDs have no prediction metrics.\nPlease double check your TMHMM file to find out if you have missing data.\n\
				The output table will be produced without SeqIDs for which metrics are absent.\n\
				You can find the offending SeqIDs in: {warnings}')
		with open(warnings, 'w') as wf:
			for e in errList:
				print(e, file=wf)


def extract_extensive(seqids: list[str], lines: list[str]) -> dict[list]:
	"""
	Takes 2 LISTs of STRINGs and an initialised DICTIONARY of LISTs.
	Populates the DICTIONARY if matches are found between the other inputs.

	Args:
		seqids (list[str]): List of IDs (strings)
		lines (list[str]): List of lines (strings) from loaded file
		dc (dict[list]): Initialised dictionary from defaultdict(list)
	"""
	dc = defaultdict(dict)
	loc = ['inside', 'outside', 'TMhelix']

	for l in lines:
		# remove whitespaces and hashmarks
		sline = l.rstrip().lstrip('# ')
		prot = sline.split("\t")[0].split(' ')[0]
		sid, p_desc = find_id(seqids, prot)

		if sid != None:
			if not silent:
				print(f'[INFO] Found prediction for {sid}')
			dc[sid].setdefault('Protein', p_desc)

			if any(val in sline for val in header_map.values()):
				for val in header_map.values():
					if val in sline:
						dc[sid][val] = sline.split(' ')[-1]

			elif any(i in sline for i in loc):
				domains = dc[sid].setdefault('Domains', [])
				for i in loc:
					if i in sline:
						dom, srt, end = re.sub(r'\s+', '-', sline).split('-')[-3:]
						domains.append(f'{dom}: {srt}-{end}')

			else:
				dc[sid]['NOTE'] = ' '.join(re.sub(r'\s+', ' ', sline).split(' ')[2:])

	for inner in dc.values():
		inner.setdefault('NOTE', 'NA') 

	return dc
'''
# --------------- these go in extract_extensive()
	# process lines of each extracted query and transform into a table
	loc = ['inside', 'outside', 'TMhelix'] # set up strings for recognition of domain line
	table = dict()
	for k,v in cleaned_file.items():
		table[k] = dict()
		lastCol = '' 
		for i in v:
			if loc[0] in i or loc[1] in i or loc[2] in i : # extract domain line
				lastCol = lastCol + i.split("\t")[-2] + ': ' + i.split("\t")[-1].lstrip().split(' ')[0] + '-' + i.split("\t")[-1].lstrip().split(' ')[-1] + '; '
				table[k]['Domains'] = lastCol
			elif ':' in i: # extract data lines
				value = i.split(' ')[-1]
				key = ' '.join(i.split(' ')[1:-1]).rstrip().rstrip(':').replace(',', '')
				table[k][key] = value
			else: # extract extra feature line (not always present)
				key = 'NOTE'
				table[k][key] = ' '.join(i.split(' ')[2:])
		if 'NOTE' not in list(table[k].keys()): # manage column for extra feature when data is absent
			table[k]['NOTE'] = 'NA'
'''


def extract_oneperline(seqids: list[str], lines: list[str], cols: dict[str]) -> dict[list]:
	"""
	Takes 2 LISTs of STRINGs and an initialised DICTIONARY of LISTs.
	Populates the DICTIONARY if matches are found between the other inputs.

	Args:
		seqids (list[str]): List of IDs (strings)
		lines (list[str]): List of lines (strings) from loaded file
		dc (dict[list]): Initialised dictionary from defaultdict(list)
	"""
	dc = dict()

	for l in lines:
		fields = l.strip().split("\t")
		sid, p_desc = find_id(seqids, fields[0])
		if sid != None:
			if not silent:
				print(f'[INFO] Found prediction for {sid}')
			dc[sid] = [ tuple(i.split('=')) if '=' in i else i for i in fields ]

	return dc


# main scripts
def run(input_path: str, idlist: str, output: str | None = None):
	"""
	Takes 3 STRINGs (path to files): reads input from STRING 1, ids from STRING 2, saves processed file to output (STRING 3).

	Args:
		input_path (str): Path to input file
		idlist (str): Path to seqids file
		output (str | None): Path to output file or None. Default: None.
	"""
	# read input file as list of lines
	with open(input_path, 'r') as fh:
		whole_file = fh.readlines()
	# read list of ids as list of strings
	with open(idlist, 'r') as ids:
		seqids = [ sid.rstrip() for sid in ids.readlines() ]

	if not whole_file:
		raise UnreadableInputError('Cannot read input file or input file is empty.')
	elif whole_file and any(line.startswith('#') for line in whole_file):
		if not silent:
			print('[INFO] \'#\' characters detected. Executing for extensive format...')
		dc = extract_extensive(seqids, whole_file)
	else:
		if not silent:
			print('[INFO] Executing for format \'One protein per line\'...')
		dc = extract_oneperline(seqids, whole_file, header_map)

	print(dc)
	sys.exit()
	if len(dc.keys()) == 0:
		raise SeqIDNotFoundError('No match for SeqIDs. Please check the input files.')

	if output:
		save_output(output, dc)

def main():
	# set default output name
	date = datetime.now()
	outName = 'tmhmm_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='tmhmm_extract.py', description='Extracts data from TMHMM plain text output and creates a tsv table with the prediction results.')
	parser.add_argument('-i', '--input', help="The plain text output copied from TMHMM's results web page - formats 'One protein per line' or 'Extensive, no images'.", required=True)
	parser.add_argument('-l', '--idlist', help='A single column, plain text table with all SeqIDs to extract.', required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output tsv table. Default: tmhmm_results_short_date_time.csv in current directory.')
	args = parser.parse_args()


	try:
		run(args.input, args.idlist, args.output)

		print(f"[INFO] Finished. Output file is {args.output}" if os.path.exists(args.output) else "[ERROR] Failed to write output file.")

	except (FileNotFoundError, ValueError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 1

'''

# ---------- these go in run()
	# check for missing data
	checkers(cleaned_file)
'''



if __name__ == "__main__":
	sys.exit(main())
