#!/usr/bin/env python

# vaxijen_extract.py
__version__ = '1.2.0'

# manage libraries
import argparse, sys, os, re
from datetime import datetime

# define functions
## saves output table in standardised format
def save_output(ofi, ret_obj):
	header = []
	with open(ofi, 'w') as out:
		for k,v in ret_obj.items():
			if len(header) == 0:
				header = list(v.keys())
				print(*header, sep="\t", file=out)
			row = [ v[i] for i in header ]
			print(*row, sep="\t", file=out)

## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj):
	save_output('vaxijen.tmp', ret_obj)


def run():
	print('WIP')


# main script
def main():
	# set default output name
	date = datetime.now()
	outName = 'vaxijen_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# manage arguments 
	parser = argparse.ArgumentParser(prog='vaxijen_extract.py', description="Extracts data from Vaxijen output (as copied from the results web page) and creates a tsv table with the prediction results.\
		Usage: python3 vaxijen_extract.py -i <input_file>")
	parser.add_argument('-i', '--input', help="Vaxijen prediction output.", required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output tsv table. Default: vaxijen_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	# read file and parse line by line to extract values for each protein
	fh = open(args.input, 'r')

	regx = r'^>[A-Za-z0-9_\.]+\|([A-Za-z0-9_]+)\|([A-Za-z0-9_]+) (.+$)' #r'^>([A-Z0-9_.]+) (.+\]) (.+[0-9]) (\(.+ \))\.$'
	dc = dict()

	for line in fh.readlines():
		if line.startswith('>'):
			stripLine = line.rstrip()
			mstring = re.match(regx, stripLine)
			seqid = mstring.group(1)
			prot, pred_v = mstring.group(3).split(' (')[0:2] #mstring.group(2)
			pred_v = pred_v.replace(')', '') #mstring.group(3)
			pred_s = re.search(r'Probable [A-Z\-]+', mstring.group(3)).group(0) #mstring.group(4)

			#dc[seqid] = {'SeqID': seqid, 'Protein': prot, pred_v.split(' = ')[0]: pred_v.split(' = ')[1], 'Prediction': pred_s.lstrip('( ').rstrip(' )')}
			dc[seqid] = {'SeqID': seqid,
						'Protein': prot,
						'Protective Antigen Prediction': pred_v,
						'Prediction': pred_s}
		else:
			continue

	fh.close()

	# kill script and throw error if input is not a valid file
	if len(dc.keys()) == 0:
		print('ERROR: could not find any row of Vaxijen output in the input file.\nPlease ensure the input file contains rows from Vaxijen\'s output as they appear on the web page of Vaxijen results.')
		sys.exit(1)

	# print table to output file
	save_output(args.output, dc)

	# print exit message and exit
	if os.path.exists(args.output):
		print("Finished. Output file is " + args.output)
	else:
		print('Failed to write output file.')

	return dc


if __name__ == "__main__":
	sys.exit(main())