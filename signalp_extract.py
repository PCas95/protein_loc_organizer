#!/usr/bin/env python

# signalp_extract.py
__version__ = '1.2.0'

# manage libraries
import argparse, json, sys, os
from datetime import datetime

# define functions
## saves output table in standardised format
def save_output(ofi, ret_obj):
	with open(ofi, 'w') as out:
		print('SeqID', 'Cleavage site', 'Likelihood SP (Sec/SPI)', 'Likelihood TAT SP (Tat/SPI)', 'Likelihood Lipoprotein SP (Sec/SPII)', 'Likelihood Other', 'Prediction', sep="\t", file=out)
		for i in ret_obj['SEQUENCES']:
			cs = ret_obj['SEQUENCES'][i]['CS_pos'] if ret_obj['SEQUENCES'][i]['CS_pos'] != '' else 'NA'
			lh = ret_obj['SEQUENCES'][i]['Likelihood']
			pr = ret_obj['SEQUENCES'][i]['Prediction']
			pt = ret_obj['SEQUENCES'][i]['Protein_types']
			print(i, cs, *lh, pr, sep="\t", file=out)

## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj):
	save_output('signalp.tmp', ret_obj)

# main script
def main():

	# set default output name
	date = datetime.now()
	outName = 'signalp5_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# manage arguments 
	parser = argparse.ArgumentParser(prog='signalp5_extract.py', description="Extracts data from SignalP5 .json output and creates a tsv table with the prediction results.\
		Usage: python3 signalp5_extract.py -i <input_json>")
	parser.add_argument('-i', '--input', help="The SignalP5 json output.", required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output tsv table. Default: signalp5_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	# load json file
	try:
		with open(args.input, 'r') as fh:
			js = json.load(fh)
	except json.decoder.JSONDecodeError:
		print('Error: input file could not be recognised as a json file.\nPlease make sure you are using the json summary file produced by SignalP as input.')
		sys.exit(1)

	# output csv table with values from desired keys
	save_output(args.output, js)

	# print exit message and exit
	if os.path.exists(args.output):
		print(f'Finished. Output file is {args.output}')
	else:
		print('Failed to write output file.')

	# return statement to pass to exp_tmp() - for when ran with Proteus wrapper
	return js


if __name__ == "__main__":
	main()
	sys.exit(0)