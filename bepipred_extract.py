#!/usr/bin/env python

# bepipred_extract.py
__version__ = '1.0.0'

# manage libraries
import argparse, os, sys, re
from datetime import datetime
from itertools import groupby
from operator import itemgetter

i = None
s = None
silent = False

def header_creator(i, s, dc):
	with open(i) as fh:
		header = fh.readline().rstrip().split(s)
		try:
			if header[0] != 'Entry' and header[1] != 'Position' and header[2] != 'AminoAcid' and header[3] != 'Exposed/Buried':
				header_error()
		except IndexError:
			header_error()

		for line in fh.readlines():
			cols = line.rstrip().split(s)
			seqid = cols[0]
			
			if seqid not in dc.keys():
				dc[seqid] = {cols[1]: { header[2]: cols[2], header[3]: cols[3],
										header[4]: cols[4], header[5]: cols[5], 
										header[6]: cols[6], header[7]: cols[7], header[8]: cols[8] }
							}
			else:
				dc[seqid][cols[1]] = { header[2]: cols[2], header[3]: cols[3],
										header[4]: cols[4], header[5]: cols[5], 
										header[6]: cols[6], header[7]: cols[7], header[8]: cols[8] }
	return header


def save_output(ofi, ret_obj, header):

	with open(ofi, 'w') as out:
		## write header to output table
		print(*header, sep=',', file=out)

		## start navigating dictionary: for each protein ID
		for k in ret_obj.keys():
			if silent == False:
				print('Processing entry ' + k + '...')

			### create filtered list of aminoacid positions exposed to solvent and with probability >=0.5
			aas = []
			for innKey in ret_obj[k].keys():
				if ret_obj[k][innKey][header[3]] == 'E' and float(ret_obj[k][innKey][header[8]]) >= 0.5:
					aas.append(int(innKey))

			### create list of consecutive aminoacids from filtered list (list of tuples)
			aas.sort()
			ranges =[]

			for z,g in groupby(enumerate(aas), lambda x: x[0] - x[1]):
				grp = (map(itemgetter(1), g))
				grp = list(map(int, grp))
				ranges.append((grp[0], grp[-1]))
			
			### extract from list of consecutive aminoacids only those that are groups of 4 or more (list of tuples)
			exposed = []
			for tup in ranges:
				if tup[-1] - tup[0] >= 3:
					exposed.append((tup[0], tup[-1]))
			if silent == False:
				print('Probable EETS found in ' + k + ': ' + str(len(exposed)))

			### for each group of 4+ consecutive aminoacids, output line with data from original file
			if len(exposed) != 0:
				for e in exposed:
					for i in range(e[0], e[-1]+1):
						vals = list(ret_obj[k][str(i)].values())
						print(k, i, *vals, sep=',', file=out)
					#print(k, exposed)
					#print([ret_obj[k][str(i)] for i in range(e[0], e[-1]+1)])


def main():

	# set default output name
	date = datetime.now()
	outName = 'bepipred_results_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.csv'
	sep = ';'

	# manage arguments 
	parser = argparse.ArgumentParser(prog='bepipred_extract.py', description="Extracts data from Bepipred output (as copied from the results web page) and creates a csv table with the prediction results.\
		Usage: python3 bepipred_extract.py -i <input_file>")
	parser.add_argument('-i', '--input', help="Bepipred prediction output.", required=True)
	parser.add_argument('-s', '--separator', default=sep, help="Field separator used in the input csv file. Default is semicolon (';', default for Bepipred output), but other separators can be specified to override this.")
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output csv table. Default: bepipred_results_date_time.csv in current directory.')
	args = parser.parse_args()

	global i,s
	i = args.input
	s = args.separator

	head_error_message = "ERROR: header of input table does not match that of Bepipred's output.\n\
	Please check that:\n\
		1. you are providing a valid Bepipred output csv table\n\
		2. your table doesn't have a poorly formatted or absent header\n\
		3. you are passing the correct field separator for this table (run bepipred_extract.py --help)"
	
	def header_error():
		print(head_error_message)
		sys.exit(1)

	# read and process file to extract aminoacid positions from each protein, the SeqIDs and the values in their columns
	dc = dict()
	innerDic = dict()

	header = header_creator(i, s, dc)

	# kill script and throw error if input is not a valid file or no SeqIDs are found
	if len(dc.keys()) == 0:
		print("ERROR: could not retrieve any sequence ID from input file.\n\
		Please ensure the input is Bepipred's csv output and that you are specifying the correct separator (default ';', check manual with " + os.path.basename(__file__) + " --help)")
		sys.exit(1)

	# extract aas exposed to solvent and create filtered output table
	# criteria: epitope probability >= 0.5 & exposed/buried == E & consecutive aminoacids >= 4
	save_output(args.output, dc, header)

	# print exit message and exit
	if os.path.exists(args.output):
		print("Finished. Output file is " + args.output)
	else:
		print('Failed to write output file.')

	return dc


def exp_tmp(ret_obj):
	global silent
	silent = True
	tmp_header = header_creator(i, s, ret_obj)
	save_output('bepipred.tmp', ret_obj, tmp_header)


if __name__ == "__main__":
	main()
	sys.exit(0)



# costruisco un dizionario come questo:		
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
