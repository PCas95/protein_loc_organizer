#!/usr/bin/env python

# list reader for SeqIDs 
def id_reader(seqID: str):
	r_seqID = seqID.replace('|', '_').split(' ')[0]
	return r_seqID