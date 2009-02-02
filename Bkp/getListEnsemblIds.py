#!/usr/bin/env python
# Legge i file contenuti in ./Ids e ottiene un oggetto con la lista di ensembl ID con una percentuale di identità maggiore di un cut-off.
# Segnala anche i casi in cui vi sono più di un ortologo per EnsemblId.


import os


class HEnsemblID:
	"""
	Human ensembl ID with its orthologs.
	"""
	pass


def readTSV(fname):
	for line in file(fname).readlines:
		if line.startswith('ENS'):
			print line



