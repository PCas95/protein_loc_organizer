#!/usr/bin/env python

# psortb_extract.py
__version__ = '1.26.01'

# libraries
import argparse, sys, os
from datetime import datetime
from common_utils import id_reader, find_id


# globals
sep = "\t"
silent = False

# functions
## saves output table in standardised format
def save_output(ofi, ret_obj):
	with open(ofi, 'a') as ofh:
		for x in ret_obj:
			print(x,file=ofh)


## parser for psortb short output
def ex_short(l: str, ids: list[str], orows: list):
	if l.startswith('SeqID'):
		cols = l.split(sep)
		if len(cols) > 3:
			raise ValueError('Number of columns does not match: malformed PSortB SHORT output.')
		new_cols = [cols[0], 'Protein', *cols[1:]]
		#orows.append("\t".join(new_cols))
		#print(*new_cols, sep="\t", file=ofh)
	else:
		fields = l.split(sep)
		seqid, desc = find_id(ids, fields[0])
		new_fields = [seqid, desc, *fields[-2:]]
		#orows.append("\t".join(new_fields))
		print(fields)
		print(new_fields[0])
		print(new_fields[1])
		print(*new_fields, sep="\t")
		#print(*new_fields, sep="\t", file=ofh)


## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj):
	save_output('psortb.tmp', ret_obj)


# main scripts
def run(input_path: str, idlist: str, output: str | None = None):
	"""
	Takes 2 STRINGs (path to files): reads input from STRING 1 and saves processed file to output (STRING 2).
	Output STRING can be None (when launched by wrapper). Default: None.
	
	Args:
		input_path (str): Path to input file
		output (str | None): Path to output file or None. Default: None.

	Returns:
		dc (dict): Extracted and cleaned data (dict object) from input file, ready for print to output table.
	"""
	
	# load seqids
	with open(idlist, 'r') as ids:
		seqids = [ sid.rstrip() for sid in ids.readlines()]
	# process table
	with open(input_path, 'r') as fh:
		orows = []
		for line in fh.readlines():
			l = line.rstrip()
			ex_short(l, seqids, orows)
	
	if output:
		save_output(args.output, orows)


def main():
	# define global variables (needed for command line arguments and processing)
	date = datetime.now()
	outName = 'psortb_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='psortb_extract.py', description="Extracts data from PSortB plain text SHORT output and creates a tsv table with the prediction results.")
	parser.add_argument('-i', '--input', help='A plain text file with PSortB SHORT output.', required=True)
	parser.add_argument('-l', '--idlist', help='A single column, plain text table with all SeqIDs to extract.', requred=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the primary output tsv table. Default: psortb_results_date_time.tsv in current directory.')
	args = parser.parse_args()


	run(args.input, args.idlist, args.output)
	# print exit message and exit
	if os.path.exists(args.output):
		print('Finished. The output table is ' + args.output)
	else:
		print('Failed to write output file ' + args.output)

	return orows


if __name__ == "__main__":
	main()
	sys.exit(0)