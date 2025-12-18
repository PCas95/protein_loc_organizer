#!/usr/bin/env python

# tmhmm_extract.py
__version__ = '1.2.0'

# manage libraries
import argparse, sys, os, re
from datetime import datetime
from collections import defaultdict

# define functions
## saves output table in standardised format
def save_output(ofi, ret_obj):
	## get full column range for output table
	nCol_max = 0
	id_max_col = ''
	for x in ret_obj.keys():
		current_nCol = len(ret_obj[x].keys())
		if current_nCol > nCol_max:
			nCol_max = current_nCol
			id_max_col = x
	## write table
	with open(ofi, 'w') as out:
		print('SeqID', ','.join(list(ret_obj[id_max_col].keys())), sep=',', file=out)
		for c in ret_obj.keys():
			values = [ret_obj[c][x] for x in list(ret_obj[id_max_col].keys())]
			print(c, *values, sep=',', file=out)

## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj):
	save_output('tmhmm.tmp', ret_obj)

## warning/error check function
def checkers(dc: dict):
	# produce warning log + message if there are no results for some SeqIDs
	errList = []
	for k in dc.keys():
		if len(dc[k]) < 6:
			errList.append(str(dc[k]))
	if len(errList) == len(dc.keys()):
		print('ERROR: It seems that none of the provided SeqIDs has TMHMM prediction metrics.\nPlease double check the input file.')
		sys.exit(1)
	elif errList == True:
		warnings = f'warnings_tmhmm_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.txt'
		print(f'WARNING: One or more SeqIDs have no prediction metrics.\nPlease double check your TMHMM file to find out if you have missing data.\n\
				The output table will be produced without SeqIDs for which metrics are absent.\n\
				You can find the offending SeqIDs in: {warnings}')
		with open(warnings, 'w') as wf:
			for e in errList:
				print(e, file=wf)


# main script
def main():
	# set default output name
	date = datetime.now()
	outName = 'tmhmm_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.csv'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='tmhmm_extract.py', description="Extracts data from TMHMM plain text output and creates a csv table with the prediction results.\
		Usage: python3 tmhmm_extract.py -i <input_file>")
	parser.add_argument('-i', '--input', help="A plain text file with TMHMM's output, as copied from the results web page.", required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output csv table. Default: tmhmm_results_short_date_time.csv in current directory.')
	args = parser.parse_args()

	# read whole file as list of lines
	with open(args.input, 'r') as fh:
		whole_file = fh.readlines()

	# set up global variables for extraction loop
	id_pattern = r'[A-Z0-9\_\.]+'
	ids = set()
	cleaned_file = defaultdict(list)

	# loop through tmhmm text file to extract queries and respective result lines
	for line in whole_file:
		sline = line.rstrip().lstrip('# ') # removes whitespaces and hashmarks
		if re.search(id_pattern, sline): # groups all lines by SeqID
			seqid = re.search(id_pattern, sline).group(0)
			ids.add(seqid)
			cleaned_file[seqid].append(sline)

	# throw error if no SeqID was found
	if len(cleaned_file.keys()) == 0:
		print('ERROR: no match for SeqIDs. Please check the input file.')
		sys.exit(1)

	# process lines of each extracted query and transform into a table
	loc = ['inside', 'outside', 'TMhelix'] # set up strings for recognition of domain line
	table = dict()
	for k,v in cleaned_file.items():
		table[k] = dict()
		lastCol = '' 
		for i in v:
			if loc[0] in i or loc[1] in i or loc[2] in i : # extract domain line
				lastCol = lastCol + i.split("\t")[-2] + ': ' + i.split("\t")[-1].lstrip().split(' ')[0] + '-' + i.split("\t")[-1].lstrip().split(' ')[-1] + '; '
				table[k]['Domains'] = lastCol
			elif ':' in i: # extract data lines
				value = i.split(' ')[-1]
				key = ' '.join(i.split(' ')[1:-1]).rstrip().rstrip(':').replace(',', '')
				table[k][key] = value
			else: # extract extra feature line (not always present)
				key = 'NOTE'
				table[k][key] = ' '.join(i.split(' ')[2:])
		if 'NOTE' not in list(table[k].keys()): # manage column for extra feature when data is absent
			table[k]['NOTE'] = 'NA'

	# check for missing data
	checkers(cleaned_file)

	# write output table
	save_output(args.output, table)

	# print exit message and exit
	if os.path.exists(args.output):
		print(f'Finished. Output file is {args.output}')
	else:
		print('Failed to write output file.')

	# return statement to pass to exp_tmp() - for when ran with Proteus wrapper
	return table


if __name__ == "__main__":
	main()
	sys.exit(0)