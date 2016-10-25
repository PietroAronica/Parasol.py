#! /usr/bin/python2.7

import os, sys
from optparse import OptionParser
from subprocess import Popen, PIPE

def run(cmd, output):
	cmdline = "tleap -f %s >> %s"%(cmd, output)
	os.system(cmdline)

def createprmtops(ff='ff99SB+', pbradii='Null', dat='Param_files/Stock/VxI_hybrid.dat', lib='VXI.lib', ion='Null', box='Null', solvbox='8.0', output='leap.log'):
	ctrl = open('lyp.in', 'w')
	ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("source %s\n"%dat)
	if pbradii != 'Null':
		ctrl.write("set default pbradii %s\n"%pbradii)
	ctrl.write("loadoff %s\n"%lib)
	ctrl.write("loadamberparams 0_100.frcmod\n")
	ctrl.write("Mutanda = loadpdb Mutabond.pdb\n")
	if ion in ['Cl-', 'Na+']:
		ctrl.write("addions Mutanda %s 0\n" %ion)
	if box != 'Null':
		ctrl.write("set Mutanda box %s\n" %box)
	if solvbox != 'Null':
		ctrl.write("solvateoct Mutanda TIP3PBOX %s\n" %solvbox)
	ctrl.write("saveamberparm Mutanda Solv_0_100.prmtop Solv_0_100.inpcrd\n")
	ctrl.write("loadamberparams 10_90.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_10_90.prmtop Solv_10_90.inpcrd\n")
	ctrl.write("loadamberparams 20_80.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_20_80.prmtop Solv_20_80.inpcrd\n")
	ctrl.write("loadamberparams 30_70.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_30_70.prmtop Solv_30_70.inpcrd\n")
	ctrl.write("loadamberparams 40_60.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_40_60.prmtop Solv_40_60.inpcrd\n")
	ctrl.write("loadamberparams 50_50.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_50_50.prmtop Solv_50_50.inpcrd\n")
	ctrl.write("loadamberparams 60_40.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_60_40.prmtop Solv_60_40.inpcrd\n")
	ctrl.write("loadamberparams 70_30.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_70_30.prmtop Solv_70_30.inpcrd\n")
	ctrl.write("loadamberparams 80_20.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_80_20.prmtop Solv_80_20.inpcrd\n")
	ctrl.write("loadamberparams 90_10.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_90_10.prmtop Solv_90_10.inpcrd\n")
	ctrl.write("loadamberparams 100_0.frcmod\n")
	ctrl.write("saveamberparm Mutanda Solv_100_0.prmtop Solv_100_0.inpcrd\n")
	ctrl.write("quit\n")
	ctrl.close()
	run('lyp.in', output)
