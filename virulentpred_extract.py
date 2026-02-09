#!/usr/bin/env python

# virulentpred_extract.py
__version__ = '2.0.0'

# manage libraries
import argparse, sys, os
import pandas as pd
from datetime import datetime

# define functions
## helper function to read csv/excel dinamically
EXCEL_EXTS = {'.xls', '.xlsx', '.ods'}

def load_table(input_path: str) -> pd.DataFrame:
	"""
	Takes a CSV/Excel file and loads it into a DataFrame.
	
	Args:
		input_path (str): path to file (csv or excel table)
	
	Returns:
		pandas.DataFrame: the table as dataframe object
	"""
	ext = input_path.split('.')[-1].lower()

	# Prefer Excel reader for known Excel extensions
	if ext in EXCEL_EXTS:
		try:
			return pd.read_excel(input_path, dtype=str)
		except Exception as e:
			raise ValueError(f"Failed to read Excel file '{input_path}': {e}")
	else:
		try:
			return pd.read_csv(input_path, sep=',', dtype=str)
		except Exception as e:
			raise ValueError(f"Failed to read csv file '{input_path}': {e}")

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
def run(input_path: str, idlist: list, output: str | None = None) -> pd.DataFrame:
	"""
	Takes 2 STRINGs (path to files): reads input from STRING 1 and saves processed file to output (STRING 2).
	Output STRING can be None (when launched by wrapper). Default: None.
	
	Args:
		input_path (str): Path to input file
		output (str | None): Path to output file or None. Default: None.

	Returns:
		dc (dict): Extracted and cleaned data (dict object) from input file, ready for print to output table.
	"""
	df = load_table(input_path)

	# check steps: throw error if input is not valid or if predictions are missing
	if df.empty or len(df.index) == 0:
		raise ValueError('Virulentpred results file is empty.')

	err_list = missing_data(df)
	if len(err_list) == len(df.index):
		raise ValueError("None of the provided SeqIDs has Virulentpred prediction metrics.")

	# standardise format 
	if df.columns[0] != 'SeqIDs':
		df.rename(columns={df.columns[0]: 'SeqID'}, inplace=True)

	### manage splitting of seqids here ###

	# write output only if asked (wrapper can pass output=None)
	if output:
		df.to_csv(output, sep=',', index=False)

	return df



def main(argv=None):
	# set default output name
	date = datetime.now()
	outName = f"virulentpred_results_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.csv"
	warnings = f"warnings_lipop_{str(date.strftime('%Y-%m-%d_%H-%M-%S'))}.txt"

	# manage arguments 
	parser = argparse.ArgumentParser(prog='virulentpred_extract.py', description="Formats Virulentpred tabular output to prepare it for join. Usage: python3 virulentpred_extract.py -i <input_txt>")
	parser.add_argument('-i', '--input', help="Virulentpred downloaded output, as csv or excel file.", required=True)
	parser.add_argument('-l', '--id_list', help="List of all SeqIDs to retrieve from file.", required=True)
	parser.add_argument('-o', '--output', default=outName, help='Optional: a file or path and file name for the output csv table. Default: virulentpred_results_date_time.csv in current directory.')
	args = parser.parse_args()

	try:
		df = run(args.input, args.id_list, args.output)
	except (FileNotFoundError, ValueError) as e:
		print(f'ERROR: {e}', file=sys.stderr)
		return 1

	# check step: warn + log if there are no results for some SeqIDs
	err_list = missing_data(df)

	if err_list:
		with open(warnings, 'w', encoding='utf-8') as wf:
			for rec in err_list:
				print(rec, file=wf)
		print(f'WARNING: Some SeqIDs have no prediction. See {warnings}')

	print(f"Finished. Output file is {args.output}" if os.path.exists(args.output) else "Failed to write output file.")
	return 0


if __name__ == "__main__":
	sys.exit(main())

