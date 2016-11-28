#!/usr/bin/python

from __future__ import division
import sys
import getopt
import datetime
import importlib
import Mutation_Modules
import Frcmod_creator
import PDBHandler
import Leapy
import Run

def makeinput(pdbfile, outfile, baa, faa, resid, atypes='Standard', newbox=8.0):
# Read PDB and load correct mutation module
	struct = PDBHandler.readpdb(pdbfile)
	Curr_mut = 'Mutation_Modules.' + baa + '_' + faa
	Curr_mut = importlib.import_module(Curr_mut)
# Determine the overall charge of the pdb to see what ion is needed
	ocharge = 0
	for res in struct.residue_list:
		if res.get_resname() in ['ASP', 'GLU', 'Cl-']:
			ocharge = ocharge - 1
		elif res.get_resname() in ['LYS', 'ARG', 'Na+']:
			ocharge = ocharge + 1
	if ocharge < 0:
		Ion = 'Na+'
	elif ocharge > 0:
		Ion = 'Cl-'
	else:
		Ion = 'Neutral'
	try:
		oldbox = struct.other_dict['Cryst1'].box_side()
		newbox = 'Null'
	except:
		oldbox = 'Null'
	if baa in ('GLU', 'ASP'):
		bal = 'PC'
	elif faa in ('GLU', 'ASP'):
		bal = 'PN'
	elif baa in ('ARG', 'LYS'):
		bal = 'PN'
	elif faa in ('ARG', 'LYS'):
		bal = 'PC'
	else:
		bal='Null'
# Make Mutated PDB (Mutabond)
	Curr_mut.makevxi(struct, 'Mutabond.pdb', resid)
# Create parameters and lib files
	Curr_mut.all_make()
	if atypes == 'Standard':
		Curr_mut.stock_add_to_all()
	else:
		Curr_mut.stock_add_to_all(**atypes)
	Curr_mut.lib_make('ff99SB+', outfile)
# Create prmtops
	Leapy.createprmtops(ff='ff99SB+', ion=Ion, box=oldbox, solvbox=newbox, output=outfile, bal=bal)
	Curr_mut.parmed_command()
	Run.make_store(oldbox)
