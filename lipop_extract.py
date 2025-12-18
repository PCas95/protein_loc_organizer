#!/usr/bin/env python

# lipop_extract.py
__version__ = '1.3.0'

# manage libraries
import argparse, sys, os, re
from datetime import datetime

# define functions
## saves output table in standardised format
def save_output(ofi: str, ret_obj: dict):
	"""
	Takes a PATH to file and a nested DICTIONARY. Prints key-values to output file as csv.
	
	Args:
		ofi (str): Path to output file
		ret_obj (dict): Nested dictionary
	"""
	with open(ofi, 'w') as out:
		print('SeqID', 'prediction', 'score', 'margin', 'cleavage', sep=',', file=out)
		for k in ret_obj.keys():
			print(ret_obj[k]['SeqID'], ret_obj[k]['prediction'], ret_obj[k]['score'], ret_obj[k]['margin'], ret_obj[k]['cleavage'], sep=',', file=out)

## function for main data extraction from file
def data_extractor(line: str, dc: dict):
	"""
	Takes a STRING and a DICTIONARY. Cleans string and splits its contents into key-values for dict.
	
	Args:
		line (str): String (line from input file)
		dc (dict): Initialised dictionary to update
	"""
	stripLine = line.rstrip().lstrip('# ') # removes whitespaces and hashmarks
	ini_vals = stripLine.split(' ') # splits line into separate values

	metrics = {'score': 'NA', 'margin': 'NA', 'cleavage': 'NA'}
	### match SeqID and extract data/prediction line
	seqstr = ini_vals[0]
	seqid = re.match(r'^([a-z]*)_([A-Z0-9]+)_([A-Z0-9]+)_', seqstr).group(2)
	pred =  ini_vals[1]
	### splits the 3 possible metrics in file
	for i in ini_vals[2:]:
		for m in metrics.keys():
			if i.startswith(m):
				metrics[m] = i.split('=')[1]
	# update dictionary with clean data from file
	metrics['SeqID'] = seqid
	metrics['prediction'] = pred
	dc[seqid] = metrics

## helper function for missing data report
def missing_data(dc: dict) -> list[str]:
	"""
	Takes a nested DICTIONARY of extracted data and outputs a list of IDs (keys) with missing data.
	
	Args:
		dc (dict): Nested dictionary of data extracted from input file. 

	Returns:
		list[str]: List of SeqIDs (as strings) with all missing data ('NA')
	"""
	err_list = [ sid for sid,rec in dc.items() if rec['score'] == rec['margin'] == rec['cleavage'] == 'NA' ]
	return err_list

## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj: dict):
	"""
	Takes a DICTIONARY and saves it running save_output() function with alternate parameters.
	
	Args:
		ret_obj (dict): Returned dictionary from main script
	"""
	save_output('lipop.tmp', ret_obj)


# main scripts
def run(input_path: str, output: str | None = None) -> dict:
	"""
	Takes 2 STRINGs (path to files): reads input from STRING 1 and saves processed file to output (STRING 2).
	Output STRING can be None (when launched by wrapper). Default: None.
	
	Args:
		input_path (str): Path to input file
		output (str | None): Path to output file or None. Default: None.

	Returns:
		dc (dict): Extracted and cleaned data (dict object) from input file, ready for print to output table.
	"""

	# dictionary initialisation + data cleaning and extraction
	dc = dict()
	with open(input_path, 'r') as fh:
		for line in fh.readlines():
			if 'Cut-off' in line or 'NO PLOT' in line or 'plot in' in line: # skips weird lines in case of alternate lipop outputs
				continue
			elif line.startswith('#'):
				data_extractor(line, dc)

	# check steps: throw error if input is not valid or if predictions are missing
	if not dc or len(dc.keys()) == 0:
		raise ValueError('Could not find LipoP results in the input file.')

	err_list = missing_data(dc)
	if len(err_list) == len(dc):
		raise ValueError("None of the provided SeqIDs has LipoP prediction metrics.")

	# write output only if asked (wrapper can pass output=None)
	if output:
		save_output(output, dc)

	return dc



def main(argv=None):
	# set default output name
	date = datetime.now()
	outName = f"lipop_results_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.csv"
	warnings = f"warnings_lipop_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.txt"

	# manage arguments 
	parser = argparse.ArgumentParser(prog='lipop_extract.py', description="Extracts data from LipoP output and creates a csv table. Usage: python3 lipop_extract.py -i <input_txt>")
	parser.add_argument('-i', '--input', help="The LipoP output, as copied from the summary web page of LipoP.", required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output csv table. Default: lipop_results_date_time.csv in current directory.')
	args = parser.parse_args()

	try:
		dc = run(args.input, args.output)
	except (FileNotFoundError, ValueError) as e:
		print(f'ERROR: {e}', file=sys.stderr)
		return 1

	# check step: warn + log if there are no results for some SeqIDs
	err_list = missing_data(dc)

	if err_list:
		with open(warnings, 'w', encoding='utf-8') as wf:
			for rec in err_list:
				print(rec, file=wf)
		print(f'WARNING: Some SeqIDs have no metrics. See {warnings}')

	print(f"Finished. Output file is {args.output}" if os.path.exists(args.output) else "Failed to write output file.")
	return 0


if __name__ == "__main__":
	sys.exit(main())