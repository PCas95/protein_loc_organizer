#!/usr/bin/env python

# psortb_extract.py
__version__ = '1.2.0'

# manage libraries
import argparse, sys, os, re
from datetime import datetime

# define functions
## saves output table in standardised format
def save_output(ofi, ret_obj):
	with open(ofi, 'a') as ofh:
		for x in ret_obj:
			print(x,file=ofh)

## exports tmp output table - to be used with Proteus wrapper
def exp_tmp(ret_obj):
	save_output('psortb.tmp', ret_obj)


# main script
def main():

	# define global variables (needed for command line arguments and processing)
	date = datetime.now()
	outName = 'psortb_results_short_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'
	#outLong = 'psortb_results_long_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'
	#tableType = 'short'
	sep = r'\s+'

	# set up arguments and help
	parser = argparse.ArgumentParser(prog='psortb_extract.py', description="Extracts data from PSortB plain text output and creates a tsv table with the prediction results.\
		Currently only PSortB's short output is accepted.\
		Usage: python3 psortb_extract.py -i <input_file>")
	parser.add_argument('-i', '--input', help="A plain text file with PSortB's output, either in long or short format (see '-m', '--mode').", required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the primary output tsv table (short). Default: psortb_results_short_date_time.tsv in current directory.')
	#parser.add_argument('-m', '--mode', default=tableType, help="Accepts either 'short' or 'long' as argument. Use to specify if the input file is PSortB's short or long table. Default: short.")
	args = parser.parse_args()

	# define functions
	## psortb short output
	def ex_short(l, orows):
		#ofh = open(args.output, 'a')
		if l.startswith('SeqID'):
			cols = re.split(sep, l)
			if len(cols) > 3:
				print("Number of columns does not match that of PSortB short output.")# Check your input and run with '--mode long' if you want to process a wider table.")
				sys.exit()
			new_cols = [cols[0], 'Protein', *cols[1:]]
			orows.append("\t".join(new_cols))
			#print(*new_cols, sep="\t", file=ofh)
		else:
			fields = re.split(sep, l)
			new_fields = [fields[0], ' '.join(fields[1:-2]), *fields[-2:]]
			orows.append("\t".join(new_fields))
			#print(*new_fields, sep="\t", file=ofh)
		#ofh.close()

	## psortb long output
	#def ex_long(l):
		#ofh = open(args.output, 'a')
		#olfh = open(outLong, 'a')
		#global col1 
		#global col2
		#global col3
		#if l.startswith('SeqID'):
			##SeqID,Protein,Localization,Score
			#cols = re.split(sep, l)
			#if len(cols) == 3:
				#print("It seems the input table is PSortB's short ouput. Please run with '--mode short' (or without the '--mode' argument) to avoid double output.\nNo output was produced for this run.")
				#sys.exit(1)
			#print(cols)
			#sys.exit(1)
			#col1 = cols.index('SeqID')
			#col2 = cols.index('Final_Localization')
			#col3 = cols.index('Final_Score')
			#new_cols_short = [cols[col1], 'Protein', 'Localization', 'Score']
			#new_cols_long = [cols[0], 'Protein', *cols[1:]]
			#print(*new_cols_short, sep=',', file=ofh)
			#print(*new_cols_long, sep=',', file=olfh)
		#else:
			#fields = re.split(sep, l)
			#new_fields_short = [fields[col1].split(' ')[0], ' '.join(fields[col1].split(' ')[1:]), fields[col2], fields[col3]]
			#new_fields_long = [fields[0].split(' ')[0], ' '.join(fields[0].split(' ')[1:]), *fields[1:]]
			#print(*new_fields_short, sep="\t", file=ofh)
			#print(*new_fields_long, sep="\t", file=olfh)
		#ofh.close()
		#olfh.close()

	# main loop to process table
	with open(args.input, 'r') as fh:
	#if args.mode == 'short':
		print('Executing for short table...')
		orows = []
		for line in fh.readlines():
			l = line.rstrip()
			ex_short(l, orows)
	#	elif args.mode == 'long':
	#		print('Executing for long table...')
	#		col1 = col2 = col3 = 0
	#		for line in fh.readlines():
	#			l = line.rstrip()
	#			ex_long(l)

	save_output(args.output, orows)
	# print exit message and exit
	#if args.mode == 'short':
	if os.path.exists(args.output):
		print('Finished. The output table is ' + args.output)
	else:
		print('Failed to write output file ' + args.output)
	#elif args.mode == 'long':
	#	if os.path.exists(args.output):
	#		print("Finished.\nThe long output table is " + outLong + "\nThe short output table is " + args.output)
	#	else:
	#		print('Failed to write output files.')

	return orows


if __name__ == "__main__":
	main()
	sys.exit(0)