#!/usr/bin/env python

# table_builder.py
__version__ = '2.26.02'

# manage libraries
import pandas as pd
import numpy as np
import argparse, sys, os.path
from datetime import datetime
from io_helpers import guess_input_files, run_tools, generate_output_name


# script and dependency versions
pandas_version = f'pandas version: {pd.__version__}'
numpy_version = f'numpy version: {np.__version__}'
argparse_version = f'argparse version: {argparse.__version__}'


def build_loc_table(cello_data: pd.DataFrame, psortb_data: pd.DataFrame,
					tmhmm_data: pd.DataFrame, signalp_data: pd.DataFrame,
					lipop_data: pd.DataFrame) -> pd.DataFrame:
	
	# cello + lipop
	df0 = cello_data[['SeqID', 'Protein', 'Prediction']]
	df0.columns = ['SeqID', 'Protein', 'CELLO']
	df1 = lipop_data[['SeqID', 'prediction']]
	df1.columns = ['SeqID', 'LIPOP']
	df0 = pd.merge(df0, df1, on='SeqID', how='inner')

	# + tmhmm
	if len(tmhmm_data.columns) > 7:
		df2 = tmhmm_data[['SeqID', 'Number of predicted TMHs']]
	else:
		df2 = tmhmm_data[['SeqID', 'PredHel']]
	
	df2.columns = ['SeqID', 'TMHMM']
	df0 = pd.merge(df0, df2, on='SeqID', how='inner')

	# + psortb
	df3 = psortb_data[['SeqID', 'Localization']]
	df3.columns = ['SeqID', 'PSORTb']
	df0 = pd.merge(df0, df3, on='SeqID', how='inner')

	# + signalp
	df4 = signalp_data[['SeqID', 'Prediction']]
	df4.columns = ['SeqID', 'SignalP-5.0']
	df0 = pd.merge(df0, df4, on='SeqID', how='inner')

	return df0


def score_assigner(df0: pd.DataFrame) -> pd.DataFrame:

	df_score = df0[['SeqID']].copy()

	# assign scores per tool
	cello_score = np.where(
    df0['CELLO'].str.contains(r'(?:^|;)Cytoplasmic(?:;|$)', regex=True, na=False), 0, 2)
	tmhmm_score = np.where(df0['TMHMM'].astype(int) >= 1, 2, 0)
	psortb_score = np.select([df0['PSORTb'] == 'Cytoplasmic',
							df0['PSORTb'] == 'Unknown'],
							[0, 1], default=2)
	lipop_score = np.where('SPII' in df0['LIPOP'], 2, 0)
	signalp_score = np.where('SPI' in df0['SignalP-5.0'], 2, 0)

	# create new score dataframe
	df_score = df_score.merge(pd.DataFrame({'SeqID': df0['SeqID'],
											'CELLO Score': cello_score,
											'PsortB Score': psortb_score,
											'TMHMM Score': tmhmm_score,
											'LIPOP Score': lipop_score,
											'SignalP Score': signalp_score}),
							on='SeqID',
							how='left')

	df_score['Final Score'] = df_score[['CELLO Score', 'PsortB Score','TMHMM Score','SignalP Score', 'LIPOP Score']].sum(axis=1)

	# add final sublocalisation column to input df
	final_sub = np.where(df_score['Final Score'].astype(int) >= 2, 'Non-Cytoplasmatic', 'Cytoplasmatic')

	df0 = df0.merge(pd.DataFrame({'SeqID': df_score['SeqID'],
								'Final Probable Sublocalization': final_sub}),
					on='SeqID',
					how='left')

	return df_score, df0


def build_ag_table(vaxijen_data: pd.DataFrame,
				   virulentpred_data: pd.DataFrame, 
				   bepipred_data: pd.DataFrame) -> pd.DataFrame:
	
	ag_table = pd.merge(vaxijen_data, virulentpred_data, bepipred_data, on='SeqID', how='inner')

	ag_table.to_csv('joined_antigen_table.tsv', sep="\t", index=False)
	# modify ag_table here

	return ag_table


# main functions
def run(cello: str | None, psortb: str | None, 
		tmhmm: str | None, signalp: str | None, 
		lipop: str | None, out_1: str, 
		vaxijen: str | None, virulentpred: str | None, 
		bepipred: str | None, out_2: str):

	# import tables from cellular localisation tools
	if cello and psortb and tmhmm and signalp and lipop:
		cello_data = pd.read_csv(cello, sep="\t", dtype=str)
		psortb_data = pd.read_csv(psortb, sep="\t", dtype=str)
		tmhmm_data = pd.read_csv(tmhmm, sep="\t", dtype=str)
		signalp_data = pd.read_csv(signalp, sep="\t", dtype=str)
		lipop_data = pd.read_csv(lipop, sep="\t", dtype=str)

		# build and write joined table
		loc_table = build_loc_table(cello_data, psortb_data, tmhmm_data, signalp_data, lipop_data)
		loc_table.to_csv('final_subloc_table.tsv', sep="\t", index=False)

		# calculate score for sublocalisation prediction, save intermediate and final tables
		score_table, loc_table = score_assigner(loc_table)
		score_table.to_csv('subloc_score_table.tsv', sep="\t", index=False)

		loc_table.to_csv(out_1, sep="\t", index=False)


	# import tables from antigen prediction tools
	if vaxijen and virulentpred and bepipred:
		vaxijen_data = pd.read_csv(vaxijen, sep="\t", dtype=str)
		virulentpred_data = pd.read_csv(virulentpred, sep="\t", dtype=str)
		bepipred_data = pd.read_csv(bepipred, sep="\t", dtype=str)

		ag_table = build_ag_table(vaxijen_data, virulentpred_data, bepipred_data)
		ag_table.to_csv(out_2, sep="\t", index=False)


def main():

	# set default output name
	date = datetime.now()
	outName1 = 'final_subloc_score_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'
	outName2 = 'final_antigen_score_' + str(date.strftime('%Y-%m-%d_%H-%M-%S')) + '.tsv'

	# manage inputs
	parser = argparse.ArgumentParser(description='table_builder2.py builds the final table of protein typing. Required inputs are tables from: cello, psortb, tmhmm, signalp, lipop, vaxijen, virulentpred, bepipred.')
	parser.add_argument('-v', '--version', action='version', version=f'table_builder2.py v{__version__}; pandas {pandas_version}; numpy {numpy_version}')
	parser.add_argument('-i', '--input-dir', help='directory containing input files; when used, do not provide individual input flags', required=False)
	parser.add_argument('-c', '--cello', help='path to cello .tsv table produced with cello_extract.py', required=False)
	parser.add_argument('-p', '--psortb', help='path to psortb .tsv table produced with psortb_extract.py', required=False)
	parser.add_argument('-t', '--tmhmm', help='path to tmhmm .tsv table produced with tmhmm_extract.py', required=False)
	parser.add_argument('-s', '--signalp', help='path to signalp5 .tsv table produced with signalp5_extract.py', required=False)
	parser.add_argument('-l', '--lipop', help='path to lipop .tsv table produced with lipop_extract.py', required=False)
	parser.add_argument('-j', '--vaxijen', help='path to vaxijen .tsv table produced with vaxijen_extract.py', required=False)
	parser.add_argument('-d', '--virulentpred', help='path to virulentpred .tsv table produced with virulentpred_extract.py', required=False)
	parser.add_argument('-b', '--bepipred', help='path to bepipred .tsv table produced with bepipred_extract.py', required=False)
	parser.add_argument('-o', '--output', help='path and name for main output table (final prediction for cellular sublocalization). Default: final_subloc_score_date_time.tsv in current directory.', default=outName1)
	parser.add_argument('-g', '--output_antigen', help='output path and name for final antigenicity table. Default: final_antigen_score_date_time.tsv in current directory.', default=outName2)

	args = parser.parse_args()
	# If an input directory is provided, disallow individual input flags and
	# populate inputs by guessing files in the directory.
	if args.input_dir:
		other_inputs = [args.cello, args.psortb, args.tmhmm, args.signalp, args.lipop, args.vaxijen, args.virulentpred, args.bepipred]
		if any(x is not None for x in other_inputs):
			parser.error('When --input-dir is used, do not provide individual input flags (cello, psortb, tmhmm, signalp, lipop, vaxijen, virulentpred, bepipred).')

		try:
			input_map = guess_input_files(args.input_dir)
		except Exception as e:
			parser.error(f'Failed to guess input files from {args.input_dir}: {e}')

		cello = input_map.get('cello')
		psortb = input_map.get('psortb')
		tmhmm = input_map.get('tmhmm')
		signalp = input_map.get('signalp')
		lipop = input_map.get('lipop')
		# Antigen-related files may not be discoverable; try bepipred if present
		bepipred = input_map.get('bepipred')
		vaxijen = input_map.get('vaxijen')
		virulentpred = input_map.get('virulentpred')
	else:
		cello = args.cello
		psortb = args.psortb
		tmhmm = args.tmhmm
		signalp = args.signalp
		lipop = args.lipop
		vaxijen = args.vaxijen
		virulentpred = args.virulentpred
		bepipred = args.bepipred

	# Print selected inputs
	print("[INFO] Selected input files:")
	print(f"  cello: {cello}")
	print(f"  psortb: {psortb}")
	print(f"  tmhmm: {tmhmm}")
	print(f"  signalp: {signalp}")
	print(f"  lipop: {lipop}")
	print(f"  vaxijen: {vaxijen}")
	print(f"  virulentpred: {virulentpred}")
	print(f"  bepipred: {bepipred}")

	run(cello, psortb, tmhmm, signalp, lipop, args.output, 
		vaxijen, virulentpred, bepipred, args.output_antigen)



if __name__ == "__main__":
	sys.exit(main())