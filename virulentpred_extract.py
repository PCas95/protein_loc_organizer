#!/usr/bin/env python

# virulentpred_extract.py
__version__ = '2.26.02'

# manage libraries
import argparse, sys, os
import pandas as pd
from datetime import datetime
from common_utils import id_reader, find_id


# classes
## custom error classes
class TableFormatError(ValueError):
	"""Malformed Virulentpred output"""

class MissingPredError(ValueError):
	"""Missing prediction metrics"""


# globals
EXCEL_EXTS = {'.xls', '.xlsx', '.ods'}
tableErrMsg = 'Dataframe empty: the provided Virulentpred output might be empty or malformed.'
sep = None

silent = False


# functions
## helper function to read csv/excel dinamically
def load_table(input_path: str) -> pd.DataFrame:
	"""
	Takes a CSV/Excel file and loads it into a DataFrame.
	
	Args:
		input_path (str): path to file (csv or excel table)
	
	Returns:
		pandas.DataFrame: the table as dataframe object
	"""
	ext = input_path.split('.')[-1].lower()

	if ext in EXCEL_EXTS:
		try:
			if not silent:
				print(f'[INFO] Reading {input_path} as Excel')
			return pd.read_excel(input_path, dtype=str)
		except Exception as e:
			raise ValueError(f"Failed to read Excel file '{input_path}': {e}")
	elif ext == 'csv':
		try:
			if not silent:
				print(f'[INFO] Reading {input_path} as CSV')
			return pd.read_csv(input_path, sep=',', dtype=str)
		except Exception as e:
			raise ValueError(f"Failed to read csv file '{input_path}': {e}")
	else:
		try:
			if not silent:
				print(f'[INFO] Reading {input_path} as TSV')
			return pd.read_csv(input_path, sep="\t", dtype=str)
		except Exception as e:
			raise ValueError(f"Failed to read file '{input_path}': {e}")


## helper function for missing data report
def missing_data(df: pd.DataFrame) -> list[str]:
	"""
	Takes a DATAFRAME and outputs a list of IDs (column 1) with missing data (column 2).
	
	Args:
		df (pandas.DataFrame): pandas dataframe object with at least 2 columns and IDs in column 1

	Returns:
		list[str]: List of SeqIDs (as strings) with all missing data ('NA' or empty string)
	"""
	# Ensure the DataFrame has at least 2 columns
	if df.shape[1] < 2:
		raise ValueError("Input DataFrame must have at least two columns.")
	col_id = df.columns[0]
	col_value = df.columns[1]
	# Select rows where the second column is NaN or empty string
	mask = df[col_value].isna() | (df[col_value].astype(str).str.strip() == "")
	# Return the IDs corresponding to missing data
	return df.loc[mask, col_id].tolist()


# main scripts
def run(input_path: str, idlist: str, sep: str | None = None, output: str | None = None) -> pd.DataFrame:
	"""
	Takes 2 STRINGs (path to files): reads input from STRING 1 and saves processed file to output (STRING 2).
	Output STRING can be None (when launched by wrapper). Default: None.
	
	Args:
		input_path (str): Path to input file
		output (str | None): Path to output file or None. Default: None.

	Returns:
		dc (dict): Extracted and cleaned data (dict object) from input file, ready for print to output table.
	"""
	# load table
	if sep is None:
		df = load_table(input_path)
	else:
		try:
			if not silent:
				print(f'[INFO] Reading {input_path} as separated by {sep}')
			return pd.read_csv(input_path, sep=sep, dtype=str)
		except Exception as e:
			raise ValueError(f"Failed to read plain text file '{input_path}' with separator {sep}: {e}")
	## remove second line ("non-tabular comment") and indices
	if df.loc[0].isna().any():
		df = df.drop(index=0)
	df.drop(df.columns[0], axis="columns", inplace=True)
	
	## normalise SeqIDs
	df['SeqID']= df.iloc[:, 0].apply(id_reader)
	df = df[['SeqID'] + [c for c in df.columns if c != 'SeqID']]

	if not silent:
		print(f'[INFO] Loaded table as pandas.DataFrame. N.of rows: {df.shape[0]}; N.of columns: {df.shape[1]}')

	# check step: error if input is not valid or if predictions are missing
	if df.empty or len(df.index) == 0:
		raise TableFormatError(tableErrMsg)

	# extract seqids based on input list
	with open(idlist, 'r') as ids:
		seqIDs = [ l.rstrip() for l in ids.readlines() ]
	
	matches = df.iloc[:, 0].apply(lambda x: find_id(seqIDs, x)[0])
	subset_df = df[matches.notna()].copy()

	if not silent:
		print(f'[INFO] Extracted table has {df.shape[0]} rows')

	# check step: error if input is not valid or if predictions are missing
	err_list = missing_data(subset_df)
	if len(err_list) == len(subset_df.index):
		raise MissingPredError('None of the provided SeqIDs has Virulentpred prediction metrics.')

	# write output if asked
	if output:
		subset_df.to_csv(output, sep="\t", index=False)

	return err_list


def main(argv=None):
	# set default output name
	date = datetime.now()
	outName = f"virulentpred_results_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.tsv"
	warnings = f"warnings_virulentpred_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.txt"

	# manage arguments 
	parser = argparse.ArgumentParser(prog='virulentpred_extract.py', description="Formats Virulentpred tabular output to prepare it for join.")
	parser.add_argument('-i', '--input', help="Virulentpred downloaded output, as csv or excel file.", required=True)
	parser.add_argument('-l', '--id_list', help="List of all SeqIDs to retrieve from file.", required=True)
	parser.add_argument('-s', '--separator', default=sep, help="String: separator used in input table. Override of default smart matching; use only for plain text table.")
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output tsv table. Default: virulentpred_results_date_time.tsv in current directory.')
	args = parser.parse_args()

	try:
		errList = run(args.input, args.id_list, args.separator, args.output)
		
		if errList:
			with open(warnings, 'w', encoding='utf-8') as wf:
				for rec in errList:
					print(rec, file=wf)
			if not silent:
				print(f'[WARNING]: Some SeqIDs have no prediction. See {warnings}')

		print(f"[INFO] Finished. Output file is {args.output}" if os.path.exists(args.output) else "[ERROR] Failed to write output file.")

	except (FileNotFoundError, ValueError) as e:
		print(f'[ERROR] {e}', file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())
