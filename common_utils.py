#!/usr/bin/env python

# list reader for SeqIDs 
def id_reader(seqID: str) -> str:
	"""
	Takes an ID with special characters as STRING and returns a normalised ID as STRING.

	Args:
		seqID (str): ID string with > or | special characters
	Returns:
		str: original string with special characters replaced ('>' -> '', '|' -> '_')
	"""
	r_seqID = seqID.split(' ')[0].replace('>', '').replace('|', '_')
	return r_seqID


# ID matcher for flattened or alternative accessions ('tr|Q8Y841|Q8Y841_LISMO' | 'tr_Q8Y841_Q8Y841_LISMO' | 'ENT65468.1')
def find_id(ids, line):
	"""
	Takes a LIST and a STRING; returns the same STRING and the item from list if the latter is exact or alternative substring.
	
	Args:
		ids (list[str]): list of ids (strings)
		line (str): line from file
	Returns:
		str: input line
		str | None: string item from input list if it's a substring of line, otherwise None
	"""
	for i in ids:

		alt = id_reader(i) if ('|' in i or i.startswith('>')) else None

		if i in line:
			return i, line
		elif alt and alt in line:
			return i, line
	
	return None, line
