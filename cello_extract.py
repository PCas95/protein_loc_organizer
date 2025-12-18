#!/usr/bin/env python

# cello_extract.py
__version__ = '1.3.0'

# manage libraries
import argparse, sys, os, re
from datetime import datetime

# define functions
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
	Takes a LIST of STRINGs and filters out non-matching strings. Returns LIST of STRINGs.
	
	Args:
		file (list): List of strings (input file read with .readlines())
		ret_obj (dict): Nested dictionary

	Returns:
		list[str]: list of filtered strings (file lines)
	"""
	for i, s in enumerate(file):
		if re.search("SeqID:", s):
			return file[i:]
	return[] # if no SeqID found, return empty list

## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj: dict):
	"""
	Takes a DICTIONARY and saves it running save_output() function with alternate parameters.
	
	Args:
		ret_obj (dict): Returned dictionary from main script
	"""
	save_output('cello.tmp', ret_obj)


# main scripts
def run(input_path: str, output: str | None = None):
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

	# process input file
	## read whole file
	with open(args.input, 'r') as fh:
		file = strip_non_data_lines(fh.readlines()) # removes any unwanted line from the top

	entries = dict() # initialise dictionary for clean lines

	current_key = None
	for line in file:
		sline = line.strip()
		if line.startswith('SeqID'):
			k = re.search(r'^SeqID: [a-zA-Z0-9\|_]*', line)
			if not k:
				continue
			current_key = k.group(0) # after removing whitespace, initialises key when SeqID is found
			entries[current_key] = []

		if not sline or sline.startswith('*****'):
			continue
		if current_key is None:
			continue
		entries[k].append(sline) # append lines after a SeqID to the that key's list

	# extract CELLO results for each query
	## prepare global variables: dictionary for output and RegExes to match lines SeqID, protein and prediction results
	tables = dict()
	## main loop: for each query result in plain text, split into lines and match only lines of interest
	for _,v in entries.items():
		for i in v:
			results = []
			### extraction of SeqID and identified protein
			if i.startswith('SeqID'):
				id_string = re.match(r'^(SeqID: ).+\|([A-Za-z0-9_]+)\|([^ ]*)', i)
				id_prot = id_string.group(2)
				tables[id_prot] = {'PROTEIN': id_prot, 'RESULTS': []}
			### extraction of prediction results
			elif re.match(r'^[A-Za-z]+\s+[0-9].+$', i):
				pred_score_loc = re.match(r'^([A-Za-z]+)\s+([0-9].+)$', i)
				loc = pred_score_loc.group(1).lstrip()
				score = pred_score_loc.group(2).split(' ')[0]

				pred = loc if pred_score_loc.group(2).split(' ')[-1] == '*' else ''
				
				tables[id_prot]['RESULTS'].append([loc, score, pred])
			### skip other lines
			else:
				continue

	if not tables or len(tables.keys()) == 0:
		raise ValueError(f"ERROR: No matched SeqIDs in input file.\nPlease check your input or run '{os.path.basename(__file__)} -h' for help.")

	if output:
		save_output(args.output, tables)

	return tables


def main():
	# set default output name
	date = datetime.now()
	outName = 'cello_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='cello_extract.py', description='Extracts data from CELLO plain text output and creates a tsv table with the prediction results.\
		Usage: python3 cello_extract.py -i <input_file> -o <output_file.csv>')
	parser.add_argument('-i', '--input', help='A plain text file with CELLO output. Can have results for multiple protein accession numbers, separated by the line of * used by CELLO.', required=True)
	parser.add_argument('-o', '--output', default=outName, help='A file or path and file name for the output csv table. Default: cello_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	try:
		tables = run(args.input, args.output):
	except (FileNotFoundError, ValueError) as e:
		print(f'ERROR: {e}', file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())
	