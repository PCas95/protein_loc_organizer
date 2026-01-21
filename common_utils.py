#!/usr/bin/env python

# list reader for SeqIDs 
def id_reader(seqID: str):
	r_seqID = seqID.split(' ')[0].replace('>', '').replace('|', '_')
	return r_seqID