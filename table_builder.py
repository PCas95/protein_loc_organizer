#!/usr/bin/env python

# table_builder.py
__version__ = '1.0.1'

# manage libraries
import pandas as pd
import numpy as np
import argparse, sys, os.path
from datetime import datetime

def main():

	# set date for output timestamp
	timestamp = str(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

	# script and dependency versions
	pandas_version = f'pandas version: {pd.__version__}'
	numpy_version = f'numpy version: {np.__version__}'
	argparse_version = f'argparse version: {argparse.__version__}'

	# manage inputs
	default_opath = './'

	parser = argparse.ArgumentParser(description='table_builder.py builds the final table of protein typing, starting from input from: cello, psortb, tmhmm, signalp, lipop, vaxijen, virulentpred, bepipred. Output tables: final_SUBLOCALIZATION and final_score_SUBLOCALIZATION')
	parser.add_argument('-v', '--version', action='version', version=f'table_builder.py v{__version__}; pandas {pandas_version}; numpy {numpy_version}')
	parser.add_argument('-c', '--cello', help='path to cello .tsv table produced with cello_extract.py', required=True)
	parser.add_argument('-p', '--psortb', help='path to psortb .tsv table produced with psortb_extract.py', required=True)
	parser.add_argument('-t', '--tmhmm', help='path to tmhmm .csv table produced with tmhmm_extract.py', required=True)
	parser.add_argument('-s', '--signalp', help='path to signalp5 .tsv table produced with signalp5_extract.py', required=True)
	parser.add_argument('-l', '--lipop', help='path to lipop .csv table produced with lipop_extract.py', required=True)
	parser.add_argument('-j', '--vaxijen', help='path to vaxijen .tsv table produced with vaxijen_extract.py', required=True)
	parser.add_argument('-d', '--virulentpred', help='path to virulentpred .tsv table produced with virulentpred_extract.py', required=True)
	parser.add_argument('-b', '--bepipred', help='path to bepipred .csv table produced with bepipred_extract.py', required=True)
	parser.add_argument('-o', '--output', help='output path. Output file will be saved at this location. By default it will be current working directory. The provided path should not contain a name for the output (pre-set)', default=default_opath)
	args = parser.parse_args()

	# import tables from localisation prediction tools
	cello_data = pd.read_csv(args.cello, sep="\t")
	psortb_data = pd.read_csv(args.psortb, sep="\t")
	tmhmm_data = pd.read_csv(args.tmhmm)
	signalp_data = pd.read_csv(args.signalp, sep="\t")
	lipop_data = pd.read_csv(args.lipop)

	# import tables from antigen prediction tools
	vaxijen_data = pd.read_csv(args.vaxijen, sep="\t")
	virulentpred_data = pd.read_csv(args.virulentpred, sep="\t")
	bepipred_data = pd.read_csv(args.bepipred)
	#mfdp_data = pd.read_csv('/home/IZSNT/p.castelli/Documents/work/Proteomics_Pipeline-17102024/mfdp2_results_2024-11-14_15-50-52.csv')

	# prepare tables for join with pd.merge
	## cello: get only rows for predicted localisation
	cello_data_ready = cello_data[cello_data['PREDICTION'] == '*']
	## rename columns
	new_cello_cols = [*cello_data_ready.columns[0:2].tolist()]

	for c in cello_data_ready.columns[2:].tolist():
		new_cello_cols.append('CELLO ' + c)

	cello_data_ready.columns = new_cello_cols

	## psortb: drop column 'Protein'
	psortb_data_ready = psortb_data.drop(columns=['Protein'])
	## rename columns
	new_psortb_cols = [psortb_data_ready.columns[0]]

	for c in psortb_data_ready.columns[1:].tolist():
		new_psortb_cols.append('Psortb ' + c)

	psortb_data_ready.columns = new_psortb_cols

	## tmhmm: keep only SeqID and N. of predicted TMHs
	cols2drop_tmhmm = ['Length', 'Exp number of AAs in TMHs', 'Exp number first 60 AAs', 'Total prob of N-in', 'Domains', 'NOTE']
	tmhmm_data_ready = tmhmm_data.drop(columns=cols2drop_tmhmm)
	## rename columns
	new_tmhmm_cols = [tmhmm_data_ready.columns[0]]

	for c in tmhmm_data_ready.columns[1:].tolist():
		new_tmhmm_cols.append('TMHMM ' + c)

	tmhmm_data_ready.columns = new_tmhmm_cols

	## signalp: keep only SeqID and Prediction
	cols2drop_signalp = ['Cleavage site', 'Likelihood SP (Sec/SPI)', 'Likelihood Other']
	signalp_data_ready = signalp_data.drop(columns=cols2drop_signalp)
	## rename columns
	new_signalp_cols = [signalp_data_ready.columns[0]]

	for c in signalp_data_ready.columns[1:].tolist():
		new_signalp_cols.append('SignalP ' + c)

	signalp_data_ready.columns = new_signalp_cols

	## lipop: keep only SeqID and Prediction
	cols2drop_lipop = ['score', 'margin', 'cleavage']
	lipop_data_ready = lipop_data.drop(columns=cols2drop_lipop)
	## rename columns
	new_lipop_cols = [lipop_data_ready.columns[0]]

	for c in lipop_data_ready.columns[1:].tolist():
		new_lipop_cols.append('Lipop ' + c)

	lipop_data_ready.columns = new_lipop_cols

	## vaxijen: keep only SeqID and Prediction
	vaxijen_data_ready = vaxijen_data.drop(columns=['Protein', 'Overall Protective Antigen Prediction'])
	## rename columns
	new_vaxijen_cols = [vaxijen_data_ready.columns[0]]

	for c in vaxijen_data_ready.columns[1:].tolist():
		new_vaxijen_cols.append('Vaxijen ' + c)

	vaxijen_data_ready.columns = new_vaxijen_cols

	## virulentpred: keep only SeqID and Prediction Results
	virulentpred_data_ready = virulentpred_data.drop(columns=['Protein Name', 'Predicted Scores'])
	## rename columns
	new_virulentpred_cols = [virulentpred_data_ready.columns[0]]

	for c in virulentpred_data_ready.columns[1:].tolist():
		new_virulentpred_cols.append('VirulentPred ' + c)

	virulentpred_data_ready.columns = new_virulentpred_cols

	## bepipred: table with 3 columns (SeqID, N. of EETS, Position of EETS)
	## initialise series for the 3 columns in new dataframe
	seqIDs = []
	n_eets = []
	epi = []
	## loop for each SeqID end get entry string + dataframe rows for that entry
	for entry, group in bepipred_data.groupby('Entry'):
		### get aa positions as a list
		positions = sorted(group['Position'].tolist())
		### extract groups of consecutive numbers from list of aa positions (returned as numpy array of arrays)
		epitopes = np.split(positions, np.where(np.diff(positions) != 1)[0]+1)
		### transform arrays into string with format 'start1-end1;start2-end2'
		epi_groups = []
		for i in epitopes:
			epi_groups.append('-'.join([str(min(i)), str(max(i))]))
		### generate full series for dataframe columns
		seqIDs.append(entry)
		n_eets.append(len(epi_groups))
		epi_string = ';'.join(epi_groups)
		epi.append(epi_string)
	## create new dataframe
	epi_dc = {'SeqID': seqIDs, 'N. of EETS': n_eets, 'Position of EETS': epi}
	bepipred_data_ready = pd.DataFrame(epi_dc, index=None)
	## rename columns
	new_bepipred_cols = [bepipred_data_ready.columns[0]]

	for c in bepipred_data_ready.columns[1:].tolist():
		new_bepipred_cols.append('BepiPred ' + c)

	bepipred_data_ready.columns = new_bepipred_cols

	# join tables
	df = pd.merge(cello_data_ready, psortb_data_ready, on='SeqID', how='inner')
	df = pd.merge(df, tmhmm_data_ready, on='SeqID', how='inner')
	df = pd.merge(df, signalp_data_ready, on='SeqID', how='inner')
	df = pd.merge(df, lipop_data_ready, on='SeqID', how='inner')
	df = pd.merge(df, vaxijen_data_ready, on='SeqID', how='inner')
	df = pd.merge(df, virulentpred_data_ready, on='SeqID', how='inner')
	df = pd.merge(df, bepipred_data_ready, on='SeqID', how='inner')

	# save table 1
	#'/home/IZSNT/p.castelli/Documents/work/Proteomics_Pipeline-17102024/final_SUBLOCALIZATION_'
#	df.to_csv(os.path.join(args.output, 'final_SUBLOCALIZATION_' + timestamp + '.csv'), sep=',', index=False)
#	df.to_excel(os.path.join(args.output, 'final_SUBLOCALIZATION_' + timestamp + '.xlsx'), index=False)

	# generate new columns for final table: scores for prediction of cytoplasmatic/membrane localisation
	## define functions for .apply()
	def cello_score(row):
		if row['CELLO LOCALIZATION'] == 'Cytoplasmic':
			return 0
		else:
			return 2

	def psortb_score(row):
		if row['Psortb Localization'] == 'Cytoplasmic':
			return 0
		elif row['Psortb Localization'] == 'Unknown':
			return 1
		else:
			return 2

	def tmhmm_score(row):
		if row['TMHMM Number of predicted TMHs'] >= 1:
			return 2
		else:
			return 0

	def signalp_score(row):
		if 'SPI' in row['SignalP Prediction']:
			return 2
		else:
			return 0

	def lipop_score(row):
		if row['Lipop prediction'] == 'SPII':
			return 2
		else:
			return 0

	def final_subl(row):
		if row['Final Score'] >= 2:
			return 'no CYT'
		else:
			return 'CYT'

	# generate new columns for final table: scores for prediction of cytoplasmatic/membrane localisation
	## create new table with final scores for subcellular localisation
	df_score = df.copy()
	df_score['CELLO Score'] = df_score.apply(cello_score, axis=1)
	df_score['PsortB Score'] = df_score.apply(psortb_score, axis=1)
	df_score['TMHMM Score'] = df_score.apply(tmhmm_score, axis=1)
	df_score['SignalP Score'] = df_score.apply(signalp_score, axis=1)
	df_score['LipoP Score'] = df_score.apply(lipop_score, axis=1)
	df_score['Final Score'] = df_score[['CELLO Score', 'PsortB Score','TMHMM Score','SignalP Score', 'LipoP Score']].sum(axis=1)
	df_score['Final Probable Sublocalization'] = df_score.apply(final_subl, axis=1)

	## additional column to table. Final probable immunogenity for cytosolic and not cytosolic proteins
	df_score['Final Probable IMMUNOGENICITY FOR CYT AND NO CYT PROTEINS'] = df_score.apply(lambda row: 'NOT IMMUNOGENIC' if row['Vaxijen Prediction'] == 'Probable NON-ANTIGEN' and row['VirulentPred Prediction results'] == 'Non-Virulent' else 'IMMUNOGENIC', axis=1)
	## additional column to table. Final probable immunogenity for cytosolic proteins only
	df_score['Final Probable IMMUNOGENICITY FOR CYT PROTEINS ONLY'] = df_score.apply(lambda row: 'IMMUNOGENIC' if row['Final Probable IMMUNOGENICITY FOR CYT AND NO CYT PROTEINS'] == 'IMMUNOGENIC' else '', axis=1)

	## drop non-score columns for subcellular localisation
	cols2drop_final = ['CELLO LOCALIZATION', 'CELLO SCORE', 'CELLO PREDICTION', 'Psortb Localization', 'Psortb Score',
						'TMHMM Number of predicted TMHs', 'SignalP Prediction', 'Lipop prediction']
	df_score.drop(columns=cols2drop_final, inplace=True)

	## reorder columns
	new_order = ['SeqID', 'PROTEIN', 'CELLO Score', 'PsortB Score', 'TMHMM Score', 'SignalP Score', 'LipoP Score', 'Final Score',
				'Final Probable Sublocalization', 'Vaxijen Prediction', 'VirulentPred Prediction results', 'BepiPred N. of EETS',
				'BepiPred Position of EETS', 'Final Probable IMMUNOGENICITY FOR CYT AND NO CYT PROTEINS', 'Final Probable IMMUNOGENICITY FOR CYT PROTEINS ONLY']
	df_score = df_score[new_order]

	# save table 2
	#/home/IZSNT/p.castelli/Documents/work/Proteomics_Pipeline-17102024/final_score_SUBLOCALIZATION_
#	df_score.to_csv(os.path.join(args.output, 'final_score_SUBLOCALIZATION_' + timestamp + '.csv'), sep=',', index=False)
#	df_score.to_excel(os.path.join(args.output, 'final_score_SUBLOCALIZATION_' + timestamp + '.xlsx'), index=False)


if __name__ == "__main__":
	main()
	sys.exit()