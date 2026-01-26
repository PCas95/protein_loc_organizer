#!/usr/bin/env python

# psortb_extract.py
__version__ = '1.26.01'

# libraries
import argparse, sys, os
from datetime import datetime
from common_utils import id_reader, find_id


# classes
## custom error classes
class PSortBFormatError(ValueError):
	"""Malformed PSortB SHORT output"""

class SeqIDNotFoundError(ValueError):
	"""No matching SeqIDs found in input"""


# globals
err_msg = f'Could not retrieve the specified SeqIDs from file.\n\
Please check the input file and the ID formats or run {os.path.basename(__file__)} -h for help.'
sep = "\t"
silent = False

# functions
## saves output table in standardised format
def save_output(ofi, ret_obj):
	with open(ofi, 'w') as ofh:
		for x in ret_obj:
			print(x,file=ofh)


## parser for psortb short output
def ex_short(l: str, ids: list[str], orows: list):
	"""
	Takes a STRING, a LIST of STRINGs and a LIST. Appends processed input to the second LIST.
	
	Args:
		l (str): line from file
		ids (list[str]): list of ids (strings)
		orows (list): empty list to which new string items will be appended
	"""
	if l.startswith('SeqID'):
		cols = l.split(sep)
		if len(cols) != 3:
			raise PSortBFormatError('Abnormal number of columns: malformed PSortB SHORT output.')
		new_cols = [cols[0], 'Protein', *cols[1:]]
		orows.append("\t".join(new_cols))
	else:
		fields = l.split(sep)
		seqid, desc = find_id(ids, fields[0])
		if seqid != None:
			new_fields = [seqid, desc, *fields[-2:]]
			orows.append("\t".join(new_fields))

			if not silent:
				print(f'[INFO] Found prediction metrics for {seqid}...')


# main scripts
def run(input_path: str, idlist: str, output: str | None = None):
	"""
	Takes 2 STRINGs (path to files): reads input from STRING 1 and saves processed file to output (STRING 2).
	Output STRING can be None (when launched by wrapper). Default: None.
	
	Args:
		input_path (str): Path to input file
		output (str | None): Path to output file or None. Default: None.
	"""
	
	# load seqids
	with open(idlist, 'r') as ids:
		seqids = [ sid.rstrip() for sid in ids.readlines()]
	# process table
	orows = []
	with open(input_path, 'r') as fh:
		for line in fh.readlines():
			l = line.rstrip()
			ex_short(l, seqids, orows)

	if len(orows) < 2:
		raise SeqIDNotFoundError(err_msg)

	if output:
		save_output(output, orows)


def main():
	# define global variables (needed for command line arguments and processing)
	date = datetime.now()
	outName = 'psortb_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='psortb_extract.py', description="Normalises PSortB TSV SHORT output and creates a tsv table with the prediction results.")
	parser.add_argument('-i', '--input', help='A plain text file with PSortB SHORT output.', required=True)
	parser.add_argument('-l', '--idlist', help='A single column, plain text table with all SeqIDs to extract.', required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the primary output tsv table. Default: psortb_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	try:
		run(args.input, args.idlist, args.output)

		print(f"[INFO] Finished. Output file is {args.output}" if os.path.exists(args.output) else "[ERROR] Failed to write output file.")

	except (FileNotFoundError, ValueError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())