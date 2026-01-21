#!/usr/bin/env python

# cello_extract.py
__version__ = '1.26.01'

# manage libraries
import argparse, sys, os, re
from datetime import datetime
from common_utils import id_reader as id_reader


# global vars
err_msg = f'[ERROR] Could not retreive the specified SeqIDs from file.\n\
Please check the input files or run {os.path.basename(__file__)} -h for help.'
silent = False


# functions

## saves output table in standardised format
def save_output(ofi: str, ret_obj: dict):
	"""
	Takes a PATH to file and a nested DICTIONARY. Prints key-values to output file as tsv.
	
	Args:
		ofi (str): Path to output file
		ret_obj (dict): Nested dictionary
	"""
	with open(ofi, 'w') as out:
			print('SeqID', 'PROTEIN', 'LOCALIZATION', 'SCORE', 'PREDICTION', sep="\t", file=out)
			for k in ret_obj.keys():
				for q,w in ret_obj[k].items():
					if q == 'RESULTS':
						for e in w:
							print(str(k), str(ret_obj[k]['PROTEIN']), str("\t".join(e)), sep="\t", file=out)


## removes empty lines / dirty lines from file
def strip_non_data_lines(file: list[str]) -> list[str]:
	"""
	Takes a LIST of STRINGs and filters out non-matching strings.
	Returns LIST of STRINGs.
	
	Args:
		file (list): List of strings (input file from .readlines())
		ret_obj (dict): Nested dictionary

	Returns:
		list[str]: list of filtered strings (file lines)
	"""
	for i, s in enumerate(file):
		if re.search("SeqID:", s):
			return file[i:]
	return [] # if no SeqID found, return empty list


## matches lines with id (exact match or alternative id)
def find_id(ids, line):
	"""
	Takes a LIST and a STRING; returns the same STRING and the item from list if the latter is exact or alternative substring.
	
	Args:
		ids (list[str]): list of ids (strings)
		line (str): line from file
	Returns:
		str: input line
		str | None: string item from input list if it's a substring of line, otherwise None
	"""
	for i in ids:
		found = False
		if i in line:
			found = True
			return i, line

		if not found:
			head, tail = i.rsplit("_", 1)
			alt = f'{head.replace("_", "|")}_{tail}'
			if alt in line:
				found = True
				return i, line
			else:
				continue
		
		if not found:
			return None, line


## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj: dict):
	# will become outdated with next wrapper release (will use run())
	"""
	Takes a DICTIONARY and saves it running save_output() function with alternate parameters.
	
	Args:
		ret_obj (dict): Returned dictionary from main script
	"""
	save_output('cello.tmp', ret_obj)


# main scripts
def run(input_path: str, idlist: str, output: str | None = None):
	"""
	Parses a CELLO plain-text output file and returns a structured dictionary.
	Takes 2 STRINGs (path to files): reads input from STRING 1 and saves processed file to output (STRING 2).
	Output STRING can be None (when launched by wrapper). Default: None.
	
	Args:
		input_path (str): Path to input file
		output (str | None): Path to output file or None. Default: None.

	Returns:
		dc (dict): Extracted and cleaned data (dict object) from input file, ready for print to output table.
	"""

	# process SeqID file
	with open(idlist, 'r') as lst:
		seqIDs = [ id_reader(l) for l in lst ]
	
	# process input file
	## read whole file
	with open(input_path, 'r') as fh:
		### remove any unwanted line from the top
		file = strip_non_data_lines(fh.readlines())

	## initialise dictionary for clean lines
	entries = dict()
	#current_key = None
	## process lines
	for line in file:
		sline = line.strip()
		### retrieve seqid line
		if sline.startswith('SeqID'):
			sline.split(' ')[1]
			k, id_string = find_id(seqIDs, sline)
			if not k:
				continue
			elif k and id_string:
				if not silent:
					print(f'[INFO] Processing entry {k}...')
				#### initialises key when SeqID is found
				entries[k] = []
				#current_key = k
		### skip empty/non data lines
		if not sline or sline.startswith('*****'):
			continue
		#if current_key is None:
			#continue
		### append data lines, separate by entry
		entries[k].append(sline) # append lines after a SeqID to the that key's list
		if not silent:
			print(f'[INFO] Found prediction results for {k}')
	# throw error if no entry matched
	if len(entries.keys()) == 0:
		raise ValueError(err_msg)


	# extract CELLO results for each query
	## prepare global variables: dictionary for output and RegExes to match lines SeqID, protein and prediction results
	tables = dict()
	## main loop: for each query result in plain text, split into lines and match only lines of interest
	if not silent:
		print(f'[INFO] Building TSV table...')
	for k,v in entries.items():
		for i in v:
			results = []
			### extraction of SeqID and identified protein
			if i.startswith('SeqID'):
				prot_desc = i.lstrip('SeqID: ')
				tables[k] = {'PROTEIN': prot_desc, 
							'RESULTS': []}
			### extraction of prediction results
			elif re.match(r'^[A-Za-z]+\s+[0-9].+$', i):
				pred_score_loc = re.match(r'^([A-Za-z]+)\s+([0-9].+)$', i)
				loc = pred_score_loc.group(1).lstrip()
				score = pred_score_loc.group(2).split(' ')[0]

				pred = loc if pred_score_loc.group(2).split(' ')[-1] == '*' else ''
				
				tables[k]['RESULTS'].append([loc, score, pred])
			### skip other lines
			else:
				continue
	if not tables or len(tables.keys()) == 0:
		raise ValueError(err_msg)

	if output:
		save_output(output, tables)

	return tables


def main():

	# set default output name
	date = datetime.now()
	outName = 'cello_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='cello_extract.py', description='Extracts data from CELLO plain text raw output and creates a TSV table with the prediction results.\
		Usage: python3 cello_extract.py -i <input_file> -o <output_file.csv>')
	parser.add_argument('-i', '--input', help='A plain text file with CELLO output. Can have results for multiple protein accession numbers, separated by the line of * used by CELLO.', required=True)
	parser.add_argument('-l', '--idlist', help='A single column, plain text table with all SeqIDs to extract.', required=True)
	parser.add_argument('-o', '--output', default=outName, help='A file or path and file name for the output csv table. Default: cello_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	try:
		tables = run(args.input, args.idlist, args.output)
		print(f"[INFO] Finished. Output file is {args.output}" if os.path.exists(args.output) else "[ERROR] Failed to write output file.")
		return 0
	except (FileNotFoundError, ValueError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())
	