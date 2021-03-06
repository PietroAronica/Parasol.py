#!/home/pietroa/.conda/envs/base/bin/python

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

def makeinput(pdbfile, outfile, baa, faa, resid, resid2='No', atypes='Standard', newbox=8.0, ff1='Standard', ff2='No', lipid='No', extralib='Null', extrafrcmod='Null', extracommand='Null', extraprep='Null', extracommand1='Null'):
# Set forcefield to be used
	if ff1 == 'Standard':
		ff1 = 'Param_files/Essentials/cmd.ff14SB+'
# If newbox not given, set to 8.0
	if newbox == 'None':
		newbox = 8.0
# Read PDB and load correct mutation module
	struct = PDBHandler.readpdb(pdbfile)
	Curr_mut = 'Mutation_Modules.' + baa + '_' + faa
	Curr_mut = importlib.import_module(Curr_mut)
# Determine the overall charge of the pdb to see what ion is needed
	ocharge = 0
	for res in struct.residue_list:
		if res.get_resname() in ['ASP', 'GLU', 'Cl-', 'ACE']:
			ocharge = ocharge - 1
		elif res.get_resname() in ['LYS', 'ARG', 'Na+', 'NME']:
			ocharge = ocharge + 1
	if ocharge < 0:
		Ion = 'Na+'
	elif ocharge > 0:
		Ion = 'Cl-'
	else:
		Ion = 'Neutral'
	try:
		oldbox = struct.other_dict['Cryst1'].box_dimensions()
		boxangle = struct.other_dict['Cryst1'].box_angle()
		newbox = 'Null'
	except:
		oldbox = 'Null'
		boxangle = 'Null'
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
	if baa == 'ARG' and faa == 'HR1':
		bal='Null'
	if baa == 'ASP' and faa == 'GLU':
		bal='Null'
	if baa == 'GLU' and (faa == 'ASP' or faa == 'HE1'):
		bal='Null'
# Determine if it is terminal
	Term = 'None'
	try:
		struct.residue_dict[resid].atom_dict['H3']
		Term = 'N'
	except:
		pass
	try:
		struct.residue_dict[resid].atom_dict['OXT']
		Term = 'C'
	except:
		pass
# Make Mutated PDB (Mutabond)
	Curr_mut.makevxi(struct, 'Mutabond.pdb', resid)
	if resid2 != 'No':
		struct = PDBHandler.readpdb('Mutabond.pdb')
		Curr_mut.makevxi(struct, 'Mutabond.pdb', resid2)
# Create parameters and lib files
	Curr_mut.all_make()
	var=Curr_mut.variablemake(sym='^')
	if Term == 'N':
			Curr_mut.stock_add_to_all_N()
	elif Term == 'C':
			Curr_mut.stock_add_to_all_C()
	else:	
		Curr_mut.stock_add_to_all(var=var)
	if Term == 'N':
		Curr_mut.lib_make_N(ff1, outfile)
	elif Term == 'C':
		Curr_mut.lib_make_C(ff1, outfile)
	else:	
		Curr_mut.lib_make(ff1, outfile, var=var)
# Create prmtops
	Leapy.createprmtops(ff1=ff1, ff2=ff2, ion=Ion, box=oldbox, solvbox=newbox, output=outfile, bal=bal, extralib=extralib, extrafrcmod=extrafrcmod, extracommand=extracommand, extraprep=extraprep, extracommand1=extracommand1)
	if Term == 'N':
		Curr_mut.parmed_command_N(lipid=lipid)
	elif Term == 'C':
		Curr_mut.parmed_command_C(lipid=lipid)
	else:	
		Curr_mut.parmed_command(lipid=lipid)
	Run.make_store(boxangle)
