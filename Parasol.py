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
                struct.other_dict['Cryst1'].box_side()
                solv = 'Yes'
        except:
                solv = 'No'
	return baa, faa, resid, solv

get_var()
baa = get_var()[0]
faa = get_var()[1]
resid = get_var()[2]
solv = get_var()[3]
pdbfile=readinput(sys.argv[1:])[1]
outfile=readinput(sys.argv[1:])[2]

# ISOLEUCINE MUTATIONS
if baa == 'ILE' and faa == 'VAL':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, solv=solv)
if baa == 'ILE' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, 'VAL', resid)
	cleanup()
	nwdir = Run.mutate_beg('/home/pietroa/Python/Store/', 'VAL', resid, solv=solv)
	os.chdir('/home/pietroa/Python')
	God.makeinput('Curr_run.pdb', outfile, 'VAL', 'ALA', resid)
	cleanup()
	Run.mutate_con('/home/pietroa/Python/Store/', 'ALA', nwdir, resid)
# VALINE MUTATIONS
if baa == 'VAL' and faa == 'ALA':
	God.makeinput(pdbfile, outfile, baa, faa, resid)
	cleanup()
	Run.mutate_beg('/home/pietroa/Python/Store/', faa, resid, solv=solv)
