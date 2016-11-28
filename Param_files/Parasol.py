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
import God
import Run
import os
import shutil

# Define input files and output files 
def readinput(argv):
	inputfile= ''
	pdbfile= ''
	outputfile= ''
	try:
		opts, args = getopt.getopt(argv,"hi:p:o:",["ifile=","pdb=","ofile="])
	except getopt.GetoptError:
		print 'God.py -i <inputfile> -p <pdbfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Usage: God.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-p", "--pdb"):
			pdbfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
	return inputfile, pdbfile, outputfile

if __name__ == "__main__":
	readinput(sys.argv[1:])

def cleanup():
        for i in range(0, 110, 10):
                os.remove('{}_{}.frcmod'.format(i, 100-i))
	os.remove('leap.log')
	os.remove('Mutabond.pdb')
	os.remove('VXI.lib')
	os.remove('lyp.in')
	os.remove('Hyb.dat')

def get_var():
	inpfile=readinput(sys.argv[1:])[0]
	pdbfile=readinput(sys.argv[1:])[1]
	outfile=readinput(sys.argv[1:])[2]
	outmake=open(outfile, 'w+')
	rightnow = str(datetime.datetime.now().strftime('DATE: %Y/%m/%d TIME: %H:%M:%S'))
	outmake.write(rightnow + '\n' )
# Define dictionary of input files
	inf = {}
# Load up input file
	f = open(inpfile)
# Parse input, assign values to variables
	data = f.readlines()
	for line in data:
	    key, value = line.split("=")
	    inf[key.strip()] = value.strip()
	f.close()
# Convert residues if necessary
	inf['residue'] = int(inf['residue'])
	resid = inf['residue']
	struct = PDBHandler.readpdb(pdbfile)
	baa = struct.residue_dict[resid].get_resname()
	faa = inf['finalaa']
	try:
		wkdir = inf['wkdir']
	except:
		wkdir = 'None'
        try:
                struct.other_dict['Cryst1'].box_side()
                solv = 'Yes'
        except:
                solv = 'No'
	return baa, faa, resid, solv, wkdir

get_var()
baa = get_var()[0]
if baa == 'HIE':
	baa = 'HIS'
faa = get_var()[1]
resid = get_var()[2]
solv = get_var()[3]
wkdir = get_var()[4]
pdbfile=readinput(sys.argv[1:])[1]
outfile=readinput(sys.argv[1:])[2]

# GLYCINE MUTATIONS
if baa == 'GLY' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'GLY' and faa in ('VAL', 'SER', 'CYS', 'THR', 'ASN', 'ASP', 'ABU'):
	God.makeinput(pdbfile, outfile, baa, 'ALA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ALA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'GLY' and faa in ('LEU', 'ILE', 'GLU', 'GLN', 'NVA', 'HCY'):
	God.makeinput(pdbfile, outfile, baa, 'ALA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ALA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'GLY' and faa in ('PRO'):
	God.makeinput(pdbfile, outfile, baa, 'ALA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ALA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', 'QUA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'QUA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'QUA', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'GLY' and faa in ('MET'):
	God.makeinput(pdbfile, outfile, baa, 'ALA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ALA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'HCY', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'HCY', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'HCY', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'GLY' and faa in ('PHE', 'TYR', 'TRP', 'HIS'):
	God.makeinput(pdbfile, outfile, baa, 'ALA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ALA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ANN', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ANN', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ANN', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'GLY' and faa in ('ARG', 'LYS'):
	God.makeinput(pdbfile, outfile, baa, 'ALA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ALA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'NLE', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NLE', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NLE', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
# ALANINE MUTATIONS
if baa == 'ALA' and faa in ('ABU', 'ASP', 'ASN', 'CYS', 'GLY', 'QUA', 'SER', 'THR', 'VAL'):
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'ALA' and faa in ('ARG', 'LYS'):
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'NLE', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NLE', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NLE', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'ALA' and faa in ('PHE', 'TYR', 'HIS', 'TRP'):
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ANN', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ANN', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ANN', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'ALA' and faa in ('ARG', 'LYS'):
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'NLE', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NLE', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NLE', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'ALA' and faa == 'PRO':
	God.makeinput(pdbfile, outfile, baa, 'QUA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'QUA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'QUA', 'PRO', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'PRO', nwdir, resid)
if baa == 'ALA' and faa in ('ILE', 'LEU', 'GLU', 'GLN'):
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'ALA' and faa in ('MET'):
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'HCY', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'HCY', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'HCY', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
# VALINE MUTATIONS
if baa == 'VAL' and faa in ('ALA', 'ASN', 'ASP', 'CYS', 'ILE', 'NVA', 'LEU', 'GLU', 'GLN', 'QUA', 'SER', 'THR'):
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'VAL' and faa == 'GLY':
	God.makeinput(pdbfile, outfile, baa, 'ALA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ALA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', 'GLY', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'GLY', nwdir, resid)
if baa == 'VAL' and faa == 'MET':
	God.makeinput(pdbfile, outfile, baa, 'NVA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'NVA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'MET', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'MET', nwdir, resid)
if baa == 'VAL' and faa == 'PRO':
	God.makeinput(pdbfile, outfile, baa, 'QUA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'QUA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'QUA', 'PRO', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'PRO', nwdir, resid)
if baa == 'VAL' and faa in ('PHE', 'HIS', 'TYR', 'TRP'):
	God.makeinput(pdbfile, outfile, baa, 'NVA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'NVA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ANN', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ANN', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ANN', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'VAL' and faa in ('ARG', 'LYS'):
	God.makeinput(pdbfile, outfile, baa, 'NVA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'NVA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'NLE', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NLE', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NLE', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
# ISOLEUCINE MUTATIONS
if baa == 'ILE' and faa in ('VAL', 'LEU', 'ABU', 'THR', 'GLU', 'GLN', 'MET'):
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'ILE' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'VAL', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'VAL', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'VAL', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
if baa == 'ILE' and faa in ('SER', 'CYS', 'ASP', 'ASN'):
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', faa, resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', faa, nwdir, resid)
if baa == 'ILE' and faa == 'GLY':
	God.makeinput(pdbfile, outfile, baa, 'VAL', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'VAL', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'VAL', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ALA', 'GLY', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'GLY', nwdir, resid)
if baa == 'ILE' and faa == 'PRO':
	God.makeinput(pdbfile, outfile, baa, 'VAL', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'VAL', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'VAL', 'QUA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'QUA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'QUA', 'PRO', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'PRO', nwdir, resid)
# METHIONINE MUTATIONS
if baa == 'MET' and faa == 'ABU':
	God.makeinput(pdbfile, outfile, baa, 'HCY', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'HCY', resid, wkdir=wkdir, solv=solv, fin='ABU')
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'HCY', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
# GLUTAMIC ACID MUTATIONS
if baa == 'GLU' and faa == 'ABU':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'GLU' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# LEUCINE MUTATIONS
if baa == 'LEU' and faa == 'ABU':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'LEU' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# THREONINE MUTATIONS
if baa == 'THR' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# GLUTAMINE MUTATIONS
if baa == 'GLN' and faa == 'ABU':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'GLN' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'ABU', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ABU', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# PROLINE MUTATIONS
if baa == 'PRO' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'QUA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'QUA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'QUA', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# TRYPTOPHAN MUTATIONS
if baa == 'TRP' and faa == 'ANN':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'TRP' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'ANN', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ANN', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ANN', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# TYROSINE MUTATIONS
if baa == 'TYR' and faa == 'ANN':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'TYR' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'ANN', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ANN', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ANN', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# QUADRINE MUTATIONS
if baa == 'QUA' and faa == 'PRO':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# PHENYLALANINE MUTATIONS
if baa == 'PHE' and faa == 'ANN':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'PHE' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'ANN', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'ANN', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ANN', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# NORVALINE MUTATIONS
if baa == 'NVA' and faa in ('ABU', 'ANN', 'NLE', 'MET'):
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# LYSINE MUTATIONS
if baa == 'LYS' and faa == 'NLE':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'LYS' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'NLE', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'NLE', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NLE', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# ARGININE MUTATIONS
if baa == 'ARG' and faa == 'NLE':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'ARG' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'NLE', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'NLE', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NLE', 'NVA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'NVA', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# NORLEUCINE MUTATIONS
if baa == 'NLE' and faa in ('ARG', 'LYS', 'NVA') :
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# NORANNULINE MUTATIONS
if baa == 'CBA' and faa == 'NVA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
if baa == 'CBA' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'NVA', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'NVA', resid, wkdir=wkdir, solv=solv, fin=faa)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'NVA', 'ABU', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ABU', nwdir, resid)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'ABU', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# ANNULINE MUTATIONS
if baa == 'ANN' and faa in ('NVA', 'PHE', 'TRP', 'TYR', 'HIS'):
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# HISTIDINE MUTATIONS
if baa == 'HIS' and faa == 'ANN':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# CYSTEINE MUTATIONS
if baa == 'CYS' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# SERINE MUTATIONS
if baa == 'SER' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# HOMOCYSTEINE MUTATIONS
if baa == 'HCY' and faa == 'MET':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# ASPARTIC ACID MUTATIONS
if baa == 'ASP' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# ASPARAGINE MUTATIONS
if baa == 'ASN' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
# AMINOBUTYRIC ACID MUTATIONS
if baa == 'ABU' and faa in ('ALA', 'CYS', 'GLU', 'GLN', 'HCY', 'ILE', 'LEU', 'NVA', 'QUA', 'SER'):
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, wkdir=wkdir, solv=solv)
