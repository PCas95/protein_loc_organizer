#!/usr/bin/env python

# signalp_extract.py
__version__ = '1.26.01'

# manage libraries
import argparse, json, sys, os
from datetime import datetime
from common_utils import id_reader, find_id


# classes
## custom error classes
class SeqIDNotFoundError(ValueError):
	"""No matching SeqIDs found in input"""


# globals
silent = False


# define functions
## saves output table in standardised format
def save_output(ofi: str, ret_obj: dict[list]):
	with open(ofi, 'w') as out:
		for k in ret_obj.keys():
			print(*ret_obj[k], sep="\t", file=out)


## parses nested json file to extract flat dictionary 
def jsn_parser(js: dict, seqids: list[str]) -> dict[list]:
	"""
	Takes a nested DICTIONARY from json import and returns a DICTIONARY of LISTs.

	Args:
		js (dict): nested dictionary from json.load()

	Returns:
		dc (dict): dictionary of lists (data lines)
	"""
	prot_types = next(iter(js["SEQUENCES"].values()))["Protein_types"]
	header = ['SeqID', 'Protein', 'Cleavage site', *prot_types, 'Prediction']

	dc = {'header': header}

	if not silent:
		print(f'[INFO] Parsing json input...')
	for i in js['SEQUENCES']:
		sid, _ = find_id(seqids, id_reader(i))
		if not silent:
			print(f'[INFO] Found prediction values for {sid}')
		cs = js['SEQUENCES'][i]['CS_pos'] if js['SEQUENCES'][i]['CS_pos'] != '' else 'NA'
		lh = js['SEQUENCES'][i]['Likelihood']
		pr = js['SEQUENCES'][i]['Prediction']
		#pt = js['SEQUENCES'][i]['Protein_types']
		dc[sid] = [sid, i, cs, *lh, pr]

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
	with open(input_path, 'r') as fh:
		js = json.load(fh)
	with open(idlist, 'r') as ids:
		seqids = [ sid.rstrip() for sid in ids.readlines() ]

	dc = jsn_parser(js, seqids)

	if len(dc.keys()) < 2:
		raise SeqIDNotFoundError('None of the listed sequence IDs was found in SignalP json file.')

	if output:
		save_output(output, dc)


def main():

	# set default output name
	date = datetime.now()
	outName = 'signalp5_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# manage arguments 
	parser = argparse.ArgumentParser(prog='signalp5_extract.py', description="Extracts data from SignalP5 .json output and creates a tsv table with the prediction results.")
	parser.add_argument('-i', '--input', help="The SignalP5 json output.", required=True)
	parser.add_argument('-l', '--idlist', help='A single column, plain text table with all SeqIDs to extract.', required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output tsv table. Default: signalp5_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	try:
		run(args.input, args.idlist, args.output)
		
		print(f"[INFO] Finished. Output file is {args.output}" if os.path.exists(args.output) else "[ERROR] Failed to write output file.")

	except json.decoder.JSONDecodeError:
		print('Malformed json file.\nPlease make sure the input is the json summary file produced by SignalP.')
		return 2
	except (FileNotFoundError, ValueError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 1	


if __name__ == "__main__":
	sys.exit(main())