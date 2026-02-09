#!/usr/bin/env python

# cello_extract.py
__version__ = '1.26.01'

# manage libraries
import argparse, sys, os, re
from datetime import datetime
from common_utils import id_reader, find_id


# classes
# custom error classes
class DataLineError(ValueError):
	"""Malformed data line with prediction metrics"""


# global vars
err_msg = f'Could not retrieve the specified SeqIDs from file.\n\
Please check the input file and the ID formats or run {os.path.basename(__file__)} -h for help.'
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
	inner_keys = set()
	for _, inner in ret_obj.items():
		if isinstance(inner, dict):
			inner_keys.update(inner.keys())

	first_col = 'SeqID'
	second_col = 'Protein'
	pred_col = 'Prediction'
	loc_cols = sorted(k for k in inner_keys if k != pred_col and k != second_col)
	cols: List[str] = [first_col, second_col] + loc_cols + [pred_col]
	
	with open(ofi, 'w') as out:
		print(*cols, sep="\t", file=out)
		for k,v in ret_obj.items():
			scores = [ v[i] for i in loc_cols if i in v.keys() ]
			output_line = [k, v[second_col], *scores, ';'.join(v[pred_col])]
			print(*output_line, sep="\t", file=out)


## removes empty/dirty lines from file
def strip_non_data_lines(file: list[str]) -> list[str]:
	"""
	Takes a LIST of STRINGs and filters out non-matching strings.
	Returns LIST of STRINGs.
	
	Args:
		file (list): List of strings (input file from .readlines())

	Returns:
		list[str]: list of filtered strings (file lines of interest)
	"""
	data = []
	for l in file:
		if l.startswith('SeqID') or re.fullmatch(r'^[A-Za-z]+\t[0-9\.\s\*]+', l):
			data.append(l)
	return data


# main scripts
def run(input_path: str, idlist: str, output: str | None = None) -> dict:

	# read from files
	with open(idlist, 'r') as lst:
		seqIDs = [ l.rstrip() for l in lst.readlines() ]

	with open(input_path, 'r') as fh:
		file = [ line.strip() for line in fh.readlines() ]
	
	# initialise variables for loop
	dc = dict()
	current_seqid = None

	# parse clean lines
	for i in strip_non_data_lines(file):
		## extract ID line
		if i.startswith('SeqID'):
			k, p_desc = find_id(seqIDs, i.lstrip('SeqID: '))
			if k is None:
				current_seqid = None
				continue

			if not silent:
				print(f'[INFO] Found data section for {k}')

			### create dictionary key-value pair and populate with protein info 
			dc[k] = dict()
			dc[k]['Protein'] = p_desc
			current_seqid = k

		else:
			## safety check
			if current_seqid is None:
				continue

			## extract data line and store in dictionary
			data_line = re.sub(r'\s+', ' ', i).split(' ')

			if len(data_line) < 2:
				raise DataLineError('Malformed prediction data line(s).')
			
			localisation, score, *extra = data_line
			pred = extra[0] if extra else None
			if not silent:
				print(f'[INFO] Processing prediction data line for {k}: {localisation}')
			dc[k][localisation] = score
			if pred is not None:
				dc[k].setdefault('Prediction', []).append(localisation)

	if output:
		save_output(output, dc)


def main():
	# set default output name
	date = datetime.now()
	outName = 'cello_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='cello_extract.py', description='Extracts data from CELLO plain text raw output and creates a TSV table with the prediction results.')
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
	