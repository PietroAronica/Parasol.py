#!/usr/bin/python

from __future__ import division
import sys
import getopt
import datetime
import importlib
import Mutation_Modules
import Frcmod_creator
import PDBHandler
import leapy

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
# Create output file and timestamp it
outfile=open(readinput(sys.argv[1:])[2], 'w+')
rightnow = str(datetime.datetime.now().strftime('DATE: %Y/%m/%d TIME: %H:%M:%S'))
outfile.write(rightnow + '\n' )
# Define dictionary of input files
inf = {}
# Load up input file
f = open(readinput(sys.argv[1:])[0])
# Parse input, assign values to variables
data = f.readlines()
for line in data:
    key, value = line.split("=")
    inf[key.strip()] = value.strip()
f.close()
# Convert residues if necessary
inf['residue'] = int(inf['residue'])
# Read PDB and load correct mutation module
struct = PDBHandler.readpdb(readinput(sys.argv[1:])[1])
beginaa = struct.residue_dict[inf['residue']].get_resname()
Curr_mut = 'Mutation_Modules.' + beginaa + '_' + inf['finalaa']
Curr_mut = importlib.import_module(Curr_mut)
try:
	Box = struct.other_dict['Cryst1'].box_side()
	inf['box'] = 'Null'
except:
	try:
		inf['box'] = float(inf['box'])
		Box = 'Null'
	except:
		inf['box'] = 8.0
		Box = 'Null'
# Make Mutated PDB (Mutabond)
Curr_mut.makevxi(struct, 'Mutabond.pdb', inf['residue'])
# Create parameters and lib files
Curr_mut.all_make()
Curr_mut.stock_add_to_all()
Curr_mut.lib_make('ff99SB+', readinput(sys.argv[1:])[2])
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
# Create prmtops
leapy.createprmtops(ff='ff99SB+', ion=Ion, box=Box, solvbox=inf['box'], output=readinput(sys.argv[1:])[2])
