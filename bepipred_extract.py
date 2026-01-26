#!/usr/bin/env python

# bepipred_extract.py
__version__ = '1.26.01'

# manage libraries
import argparse, os, sys, re
from datetime import datetime
from itertools import groupby
from operator import itemgetter


# classes
## custom error classes
class HeaderFormatError(ValueError):
	"""Malformed PSortB SHORT output"""

class SeqIDNotFoundError(ValueError):
	"""No matching SeqIDs found in input"""


# global vars
i = None
s = None
silent = False

sep = ','
outSep = "\t"
description = "bepipred_extract.py extracts epitopes exposed to solvent from Bepipred output and creates a TSV table with the prediction results.\nThe epitope is considered exposed if: 1) epitope probability >= 0.5; 2) exposed/buried == E; 3) >= 4 consecutive aminoacids meet conditions 1 and 2."
head_err_msg = f'Malformed Bepipred header. Please check your input or run {os.path.basename(__file__)} -h for help.'
seq_err_msg = f"Could not retrieve any sequence ID from input file.\n\
	Please double check input (Bepipred's csv output) and separator (default ',') or run {os.path.basename(__file__)} -h for help."


# functions
def initialiser(i: str, s: str, dc: dict) -> list:
	"""
	Takes a PATH to file, a STRING and an empty DICTIONARY. Returns a header as LIST and updates the DICTIONARY.
	
	Args:
		i (str): Path to input file
		s (str): Field separator for input file
		dc (dict): Dictionary
	Returns:
		list[str]: List of column names
	"""
	with open(i) as fh:
		header = fh.readline().rstrip().split(s)
		try:
			if header[0] != 'Entry' or header[1] != 'Position' or header[2] != 'AminoAcid' or header[3] != 'Exposed/Buried':
				raise HeaderFormatError(head_err_msg)
		except IndexError:
			raise HeaderFormatError(head_err_msg)

		for line in fh.readlines():
			cols = line.rstrip().split(s)
			seqid = cols[0]
			
			if seqid not in dc.keys():
				dc[seqid] = {cols[1]: { header[2]: cols[2], header[3]: cols[3],
										header[4]: cols[4], header[5]: cols[5], 
										header[6]: cols[6], header[7]: cols[7],
										header[8]: cols[8] }
							}
			else:
				dc[seqid][cols[1]] =  { header[2]: cols[2], header[3]: cols[3],
										header[4]: cols[4], header[5]: cols[5], 
										header[6]: cols[6], header[7]: cols[7],
										header[8]: cols[8] }
	return header


def exposed_parser(ret_obj: dict, header: list, k: str) -> list:
	"""
	Takes a nested DICTIONARY, a LIST and a STRING.
	Filters the dictionary and returns a LIST of consecutive positions (TUPLE of INTEGERS).

	Args:
		ret_obj (dict): Nested dictionary
		header (list): Header of original file to access inner keys of ret_obj
		k (str): Dictionary key
	"""

	# create filtered list of aminoacid positions exposed to solvent and with probability >=0.5
	aas = []
	for innKey in ret_obj[k].keys():
		if ret_obj[k][innKey][header[3]] == 'E' and float(ret_obj[k][innKey][header[8]]) >= 0.5:
			aas.append(int(innKey))

	# create list of consecutive aminoacids from filtered list (list of tuples)
	aas.sort()
	ranges =[]

	for z,g in groupby(enumerate(aas), lambda x: x[0] - x[1]):
		grp = (map(itemgetter(1), g))
		grp = list(map(int, grp))
		ranges.append((grp[0], grp[-1]))
	
	# extract from list of consecutive aminoacids only those that are groups of 4 or more (list of tuples)
	exposed = []
	for tup in ranges:
		if tup[-1] - tup[0] >= 3:
			exposed.append((tup[0], tup[-1]))

	return exposed


def save_output(ofi: str, ret_obj: dict, header: list):
	"""
	Takes a PATH to output file, a nested DICTIONARY and a LIST.
	Calls dictionary filter function and prints lines to file at PATH.
	
	Args:
		ofi (str): Path to output file
		ret_obj (dict): Nested dictionary (read)
		header (list[str]): List of column names
	"""
	n_lines = 0

	with open(ofi, 'w') as out:
		## write header to output table
		print(*header, sep="\t", file=out)
		
		## start navigating dictionary: for each protein ID
		for k in ret_obj.keys():
			if silent == False:
				print('[INFO] Processing entry ' + k + '...')

			exposed = exposed_parser(ret_obj, header, k)

			if silent == False:
				print('[INFO] Probable EETS found in ' + k + ': ' + str(len(exposed)))

			# for each group of 4+ consecutive aminoacids, output line with data from original file
			if len(exposed) != 0:
				for e in exposed:
					for i in range(e[0], e[-1]+1):
						vals = list(ret_obj[k][str(i)].values())
						print(k, i, *vals, sep="\t", file=out)
						n_lines += 1

	return n_lines


# main functions
def run(input_path: str, sep: str, output: str | None = None) -> dict:

	# read and process file to extract aminoacid positions from each protein, the SeqIDs and the values in their columns
	dc = dict()
	innerDic = dict()

	header = initialiser(input_path, sep, dc)

	# kill script and throw error if input is not a valid file or no SeqIDs are found
	if len(dc.keys()) == 0:
		raise SeqIDNotFoundError(seq_err_msg)

	# extract AAs exposed to solvent and create filtered output table
	# criteria: epitope probability >= 0.5 & exposed/buried == E & consecutive aminoacids >= 4
	print(f"\n{description}\n")

	# write output only if asked (wrapper can pass output=None)
	if output:
		n_lines = save_output(output, dc, header)

	if n_lines < 1:
		print("\n[WARN] No epitope exposed to solvent was found\n")

	return dc


def main():

	# set default vars
	date = datetime.now()
	outName = 'bepipred_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# manage arguments 
	parser = argparse.ArgumentParser(prog='bepipred_extract.py', description=description)
	parser.add_argument('-i', '--input', help="Raw prediction output from Bepipred.", required=True)
	parser.add_argument('-s', '--separator', default=sep, help="Field separator used in the input file. Default is comma (',', default for Bepipred output).")
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output csv table. Default: bepipred_results_date_time.csv in current directory.')
	args = parser.parse_args()

	try:
		run(args.input, args.separator, args.output)
		print(f"[INFO] Finished. Output file is {args.output}" if os.path.exists(args.output) else "[ERROR] Failed to write output file.")
		return 0
	except (FileNotFoundError, ValueError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())




# building a nested dictionary like this:
	#{
	#	XX_XXXX_X: {
	#				pos1:
	#					{
	#						exp: y/n,
	#						prob: >|<0.5
	#						...
	#					},
	#				pos2:
	#					{
	#						exp: y/n,
	#						prob: >|<0.5
	#						...
	#					}
	#				...
	#	},
	#	YY_YYYY_Y: {
	#				pos1:
	#					{
	#						exp: y/n,
	#						prob: >|<0.5
	#						...
	#					},
	#				pos2:
	#					{
	#						exp: y/n,
	#						prob: >|<0.5
	#						...
	#					}
	#				...
	#	}
	#	...
	#}
