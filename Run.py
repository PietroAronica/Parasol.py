#!/usr/bin/python

import os, sys
import shutil
import datetime
import pytraj as pt
import PDBHandler
simpath='/home/pietroa/Python/SimFiles/'

def make_store(boxangle):
	try:
		os.mkdir('Store')
	except:
		pass
	for i in range(0, 110, 10):
		shutil.move('Solv_{}_{}.prmtop'.format(i, 100-i), 'Store')
	shutil.move('Solv_0_100.inpcrd', 'Store')
	for i in range(10, 110, 10):
		os.remove('Solv_{}_{}.inpcrd'.format(i, 100-i))
	if boxangle != 90.0:
		os.chdir('Store')
		os.system('sed -i \'$s/\ 90.0000000/109.4712190/g\' Solv_*.inpcrd')
		os.system('sed -i s/9.00000000E+01/1.09471219E+02/g Solv_*.prmtop')
		os.chdir('..')

def mutate_beg(store, faa, resnum, dct='/home/pietroa/Work', wkdir='None', solv = 'No', fin=None, strip = 'No', lipid = 'No'):
	os.chdir(dct)
	if wkdir == 'None':
		if fin is None:
			fin = faa
		wkdir=datetime.datetime.now().strftime('Mutation_%m-%d_%H:%M_to{}{}'.format(fin, resnum))
	try:
		os.mkdir(wkdir)
	except:
		pass
	os.chdir(wkdir)
	nwdir = os.getcwd()
	os.mkdir('To{}'.format(faa))
	os.chdir('To{}'.format(faa))
	shutil.move(store, os.getcwd())
	if solv == 'No':
		Input = simpath + 'Method.x'
	elif solv == 'Yes':
		Input = simpath + 'Method_proceed.x'
	if lipid == 'Yes':
		Input = simpath + 'Method_lipid_proceed.x'
	os.system(Input)
	os.chdir('Run_1')
	traj = pt.load('100_0.restart', top='../Store/Solv_0_100.prmtop')
	if strip != 'No':
		traj.strip(':WAT,Cl-,Na+')
	traj.autoimage()
	pt.write_traj('Run_1.pdb', traj, overwrite=True)
	s = PDBHandler.readpdb('Run_1.pdb')
	PDBHandler.chop(s, faa, resnum)
	shutil.copyfile('Mut_leap.pdb', '/home/pietroa/Python/Curr_run.pdb')
	return nwdir

def mutate_con(store, faa, wkdir, resnum, strip = 'No', lipid = 'No'):
	os.chdir(wkdir)
	os.mkdir('To{}'.format(faa))
	os.chdir('To{}'.format(faa))
	shutil.move(store, os.getcwd())
	Input = simpath + 'Method_proceed.x'
	if lipid == 'Yes':
		Input = simpath + 'Method_lipid_proceed.x'
	os.system(Input)
	os.chdir('Run_1')
	traj = pt.load('100_0.restart', top='../Store/Solv_0_100.prmtop')
	if strip != 'No':
		traj.strip(':WAT,Cl-,Na+')
	traj.autoimage()
	pt.write_traj('Run_1.pdb', traj, overwrite=True)
	s = PDBHandler.readpdb('Run_1.pdb')
	PDBHandler.chop(s, faa, resnum)
	shutil.copyfile('Mut_leap.pdb', '/home/pietroa/Python/Curr_run.pdb')

def staple(store, resnum, dct='/home/pietroa/Work', wkdir='None'):
	os.chdir(dct)
	if wkdir == 'None':
		wkdir=datetime.datetime.now().strftime('Staple_%m-%d_%H:%M_{}-{}'.format(resnum, resnum + 4))
	try:
		os.mkdir(wkdir)
	except:
		pass
	os.chdir(wkdir)
	shutil.move(store, os.getcwd())
	Input = simpath + 'Method_proceed.x'
	os.system(Input)
	os.chdir('Run_1')
	traj = pt.load('100_0.restart', top='../Store/Solv_0_100.prmtop')
	traj.autoimage()
	pt.write_traj('Run_1.pdb', traj, overwrite=True)
	s = PDBHandler.readpdb('Run_1.pdb')
	PDBHandler.chop(s, 'NLS', resnum)
	PDBHandler.chop(s, 'NLS', resnum + 4)
