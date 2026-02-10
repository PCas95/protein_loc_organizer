#!/usr/bin/env python

# vaxijen_extract.py
__version__ = '1.26.02'

# manage libraries
import argparse, sys, os, re
from datetime import datetime
from common_utils import id_reader, find_id


# classes
## custom error classes
class MissingPredError(ValueError):
	"""Missing prediction metrics"""


# globals
silent = False


# functions
## saves output table in standardised format
def save_output(ofi: str, ret_obj: dict[dict], header: list[str]):
	with open(ofi, 'w') as out:
		print(*header, sep="\t", file=out)
		for _,v in ret_obj.items():
			row = [ v[i] for i in header ]
			print(*row, sep="\t", file=out)

## removes non-data lines
def strip_empty_lines(file: list[str]) -> list[str]:
	
	clean_file = [ line.rstrip() for line in file if line.strip() ]	
	return clean_file


# main functions
def run(input_path: str, idlist: str, output: str | None = None):

	with open(input_path, 'r') as fh:
		clean_file = strip_empty_lines(fh.readlines())
	
	if not clean_file:
		raise MissingPredError(f'File {input_path} is empty or corrupted.')

	regx = r'^(.+)(Overall Protective Antigen Prediction) =(.+)$'
	score_header = 'Overall Protective Antigen Prediction'
	pred_header = 'Prediction'
	seqid_header = 'SeqID'
	p_desc_header = 'Protein'
	header = [seqid_header, p_desc_header, score_header, pred_header]

	dc = dict()
	
	for i in clean_file:
		groups = re.fullmatch(regx, i)

		p_desc = groups[1].strip()
		seqid = id_reader(p_desc)
		g3 = groups[3].strip()
		score = g3.split(' ')[0]
		pred = re.search(r'Probable [A-Z\-]+', g3).group()

		dc.setdefault(seqid,
					{'SeqID': seqid,
					 'Protein': p_desc,
					 'Overall Protective Antigen Prediction': score,
					 'Prediction': pred}
					)

	if not silent:
		print(f'[INFO] Found {len(dc.keys())} prediction lines')

	if len(dc.keys()) == 0:
		raise MissingPredError('Could not find any prediction lines')

	with open(idlist, 'r') as ids:
		seqIDs = [ l.rstrip() for l in ids.readlines() ]

	# extract list-wise
	subset = { id_reader(i) for i in seqIDs }
	dc = { k: v for k, v in dc.items() if k in subset }

	if not silent:
		print(f'[INFO] Extracted {len(dc.keys())} matching proteins')

	if len(dc.keys()) == 0:
		raise MissingPredError('No prediction lines matching IDs')

	# write output if asked
	if output:
		save_output(output, dc, header)

	return dc


# main script
def main():
	# set default output name
	date = datetime.now()
	outName = 'vaxijen_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# manage arguments 
	parser = argparse.ArgumentParser(prog='vaxijen_extract.py', description="Extracts data from Vaxijen output (as copied from the results web page) and creates a tsv table with the prediction results.")
	parser.add_argument('-i', '--input', help="Vaxijen prediction output.", required=True)
	parser.add_argument('-l', '--id_list', help="List of all SeqIDs to retrieve from file.", required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output tsv table. Default: vaxijen_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	try:
		run(args.input, args.id_list, args.output)
	except (FileNotFoundError, ValueError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 1

	# print exit message and exit
	print(f"[INFO] Finished. Output file is {args.output}" if os.path.exists(args.output) else "[ERROR] Failed to write output file.")


if __name__ == "__main__":
	sys.exit(main())
