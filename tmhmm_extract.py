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

class ValNotFoundError(ValueError):
	"""No metric found for SeqIDs"""


# globals
silent = False
header_map = {'len': 'Length',
			'PredHel': 'Number of predicted TMHs',
			'ExpAA': 'Exp number of AAs in TMHs',
			'First60': 'Exp number, first 60 AAs'}


# functions
## saves output table in standardised format
def save_output(ofi: str, ret_obj: list[list]):
	
	with open(ofi, 'w') as out:
		for i in ret_obj:
			row = [ ';'.join(x) if isinstance(x, list) else x for x in i ]
			print(*row, sep="\t", file=out)


def checkers(lst: list[list]):
	"""
	Takes a LIST of LISTs, checks contents and raises custom errors if data is missing.

	Args:
		ls (list[list]): List of lists (rows, values in columns)  
	"""

	if len(lst) <= 2:
		raise SeqIDNotFoundError('No match for SeqIDs. Please check the input files.')

	errList = [ i[0] for i in lst if len(i) < 8 ]

	if len(errList) == len(lst) - 1:
		raise ValNotFoundError('TMHMM prediction metrics could not be found for the provided SeqIDs. Please double check the input file.')

	elif errList == True:
		warnings = f'warnings_tmhmm_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.txt'

		if not silent:
			print(f'[WARNING] Some SeqIDs have no prediction metrics.\nPlease check input file for missing data.\n\
				SeqIDs without metrics will not be present in the output table.\n\
				Offending SeqIDs: {warnings}')

		with open(warnings, 'w') as wf:
			for e in errList:
				print(e, file=wf)


def d_processor(token: str, dc: dict) -> list[list]:
	"""
	Takes a STRING and a nested DICTIONARY. Processes DICTIONARY depending on string token and returns a LIST of lines.

	Args:
		token (str): 'one' -> from oneperline format; 'ext' -> from extensive format
		dc (dict): dict[dict] from extensive or dict[list] 
		output (str | None): Path to output file or None. Default: None.
	"""
	o_lines = []
	if token == 'ext':

		tmp_inn = next(iter(dc.values()))
		header = ['SeqID'] + [ item for item in tmp_inn.keys() if isinstance(tmp_inn, dict) ]
		o_lines.append(header)

		for k,v in dc.items():
			line = [k] + [ v[i] for i in header[1:] ]
			o_lines.append(line)

	elif token == 'one':
		
		header = ['SeqID', 'Protein'] + [ key for item in next(iter(dc.values())) if isinstance(item, tuple) for key in [item[0]] ]
		o_lines.append(header)
		
		for k,v in dc.items():
			line = [k, v[0]] + [ i[-1] for i in v if isinstance(i, tuple) ]
			o_lines.append(line)

	return o_lines


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
				notes = dc[sid].setdefault('NOTE', [])
				extracted = ' '.join(re.sub(r'\s+', ' ', sline).split(' ')[2:])
				if extracted:
					notes.append(extracted)

	for inner in dc.values():
		if not inner.get('NOTE'):
			inner['NOTE'] = ['NA']

	return 'ext', dc


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

	return 'one', dc


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
		token, dc = extract_extensive(seqids, whole_file)
	else:
		if not silent:
			print('[INFO] Executing for format \'One protein per line\'...')
		token, dc = extract_oneperline(seqids, whole_file, header_map)

	out_table = d_processor(token, dc)
	checkers(out_table)

	if output:
		save_output(output, out_table)


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
	except (IndexError, KeyError, TypeError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 2


if __name__ == "__main__":
	sys.exit(main())
