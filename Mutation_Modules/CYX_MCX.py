# CYX to MCX Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/CYX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/MCX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		i = float(i)
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@SG'.format(vxi), fc['SG']/10*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), bc['SG']+((fc['CD']-bc['SG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), fc['HD2']/10*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), fc['HD3']/10*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		#print printDetails(parm, ':VXI')
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_N(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/NCYX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/NMCX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		i = float(i)
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H1'.format(vxi), bc['H1']+((fc['H1']-bc['H1'])/10)*i).execute()
                change(parm, 'charge', ':{}@H2'.format(vxi), bc['H2']+((fc['H2']-bc['H2'])/10)*i).execute()
                change(parm, 'charge', ':{}@H3'.format(vxi), bc['H3']+((fc['H3']-bc['H3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@SG'.format(vxi), fc['SG']/10*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), bc['SG']+((fc['CD']-bc['SG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), fc['HD2']/10*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), fc['HD3']/10*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		#print printDetails(parm, ':VXI')
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_C(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/CCYX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/CMCX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		i = float(i)
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@SG'.format(vxi), fc['SG']/10*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), bc['SG']+((fc['CD']-bc['SG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), fc['HD2']/10*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), fc['HD3']/10*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
                change(parm, 'charge', ':{}@OXT'.format(vxi), bc['OXT']+((fc['OXT']-bc['OXT'])/10)*i).execute()
		#print printDetails(parm, ':VXI')
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CB = struct.residue_dict[aa].atom_dict['CB']
        SG = struct.residue_dict[aa].atom_dict['SG']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HB3' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('SG', CB, SG))
                        	pdb.write(atom.halfway_between1('HD2', CB, SG))
                        	pdb.write(atom.halfway_between2('HD3', CB, SG))
			elif atom.get_name() == 'SG' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CD'))
			else:
                        	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def variablemake(sym='^'):
	var1 = sym + '1'
	var2 = sym + '2'
	var3 = sym + '3'
	var4 = sym + '4'
	var5 = sym + '5'
	var6 = sym + '6'
	var7 = sym + '7'
	var8 = sym + '8'
	var9 = sym + '9'
	var10 = sym + '0'
	var11 = sym + 'a'
	var12 = sym + 'b'
	var13 = sym + 'c'
	var14 = sym + 'd'
	var15 = sym + 'e'
	return var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, var15

def lib_make(ff, outputfile, vxi='VXI', var=variablemake()):
	intsul = var[0]
	newhyd = var[1]
	carsul = var[2]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/CYX-MCX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "S"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "SG"\n'%vxi)
	ctrl.write('set %s.1.9 name "CD"\n'%vxi)
	ctrl.write('set %s.1.10 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.12 name "C"\n'%vxi)
	ctrl.write('set %s.1.13 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "H1"\n'%vxi)
	ctrl.write('set %s.1.7 type "H1"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, intsul))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, carsul))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.12 type "C"\n'%vxi)
	ctrl.write('set %s.1.13 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect2 %s.1.9\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_N(ff, outputfile, vxi='VXI', var=variablemake()):
	intsul = var[0]
	newhyd = var[1]
	carsul = var[2]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/NCYX-NMCX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "H"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "C"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "S"\n'%vxi)
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H1"\n'%vxi)
	ctrl.write('set %s.1.3 name "H2"\n'%vxi)
	ctrl.write('set %s.1.4 name "H3"\n'%vxi)
	ctrl.write('set %s.1.5 name "CA"\n'%vxi)
	ctrl.write('set %s.1.6 name "HA"\n'%vxi)
	ctrl.write('set %s.1.7 name "CB"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.9 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.10 name "SG"\n'%vxi)
	ctrl.write('set %s.1.11 name "CD"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.13 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.14 name "C"\n'%vxi)
	ctrl.write('set %s.1.15 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N3"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "H"\n'%vxi)
	ctrl.write('set %s.1.4 type "H"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HP"\n'%vxi)
	ctrl.write('set %s.1.7 type "CT"\n'%vxi)
	ctrl.write('set %s.1.8 type "H1"\n'%vxi)
	ctrl.write('set %s.1.9 type "H1"\n'%vxi)
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, intsul))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, carsul))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.14 type "C"\n'%vxi)
	ctrl.write('set %s.1.15 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.11\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_C(ff, outputfile, vxi='VXI', var=variablemake()):
	intsul = var[0]
	newhyd = var[1]
	carsul = var[2]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/CCYX-CMCX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "S"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "O"\n'%vxi)
	ctrl.write('set %s.1.14 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "SG"\n'%vxi)
	ctrl.write('set %s.1.9 name "CD"\n'%vxi)
	ctrl.write('set %s.1.10 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.12 name "C"\n'%vxi)
	ctrl.write('set %s.1.13 name "O"\n'%vxi)
	ctrl.write('set %s.1.14 name "OXT"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "H1"\n'%vxi)
	ctrl.write('set %s.1.7 type "H1"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, intsul))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, carsul))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.12 type "C"\n'%vxi)
	ctrl.write('set %s.1.13 type "O2"\n'%vxi)
	ctrl.write('set %s.1.14 type "O2"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect2 %s.1.9\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def all_make():
	for i in range(0,110,10):
		Frcmod_creator.make ('{}_{}.frcmod'.format(i, 100-i))

def cal(x, y, i):
	num = x+((y-x)/10)*i
	return num

def lac(y, x, i):
	num = x+((y-x)/10)*i
	return num

def stock_add_to_all(var=variablemake()):
	intsul = var[0]
	newhyd = var[1]
	carsul = var[2]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(intsul, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(newhyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(carsul, 'S', 'sp3')
	p = {}
	with open('Param_files/Stock/Stock.param', 'r') as b:
		data = b.readlines()[1:]
	for line in data:
		p[line.split()[0]] = []
		for point in line.split()[1:]:
			p[line.split()[0]].append(float(point))
	b.close()
	for i in range(11):
		a = i*10
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), intsul, cal(p['0_S'][0], p['SH'][0], i), cal(p['0_S'][1], p['SH'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), carsul, cal(p['SH'][0], p['CT'][0], i), cal(p['SH'][1], p['CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', intsul), cal(p['CT_mS'][0], p['CT_S'][0], i), cal(p['CT_mS'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(intsul, carsul), cal(p['CT_mS'][0], p['CT_S'][0], i), cal(p['CT_mS'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carsul, newhyd), cal(p['CT_mSH'][0], p['CT_HC'][0], i), cal(p['CT_mSH'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carsul, 'S '), cal(p['S_S'][0], p['CT_S'][0], i), cal(p['S_S'][1], p['CT_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', intsul), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', intsul), cal(p['C_C_S'][0], p['C_C_S'][0], i), cal(p['C_C_S'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intsul, carsul), cal(p['Dritt'][0], p['C_S_C'][0], i), cal(p['Dritt'][1], p['C_S_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intsul, carsul, newhyd), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, carsul, newhyd), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intsul, carsul, 'S '), cal(p['C_S_S'][0], p['S_C_S'][0], i), cal(p['C_S_S'][1], p['S_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, carsul, 'S '), cal(p['C_S_S'][0], p['C_C_H'][0], i), cal(p['C_S_S'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carsul, 'S ', 'CT'), cal(p['C_S_S'][0], p['C_S_S'][0], i), cal(p['C_S_S'][1], p['C_S_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carsul, 'S ', '2C'), cal(p['C_S_S'][0], p['C_S_S'][0], i), cal(p['C_S_S'][1], p['C_S_S'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'S ', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'S ', '2C'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'S ', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'S ', '2C'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', intsul, carsul), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', intsul, carsul), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intsul, carsul, 'S '), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intsul, carsul, newhyd), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carsul, 'S ', 'S ', carsul), cal(p['0_1'][0], p['C_S_S_C_2'][0], i), cal(p['0_1'][1], p['C_S_S_C_2'][1], i), cal(p['0_1'][2], p['C_S_S_C_2'][2], i), cal(p['0_1'][3], p['C_S_S_C_2'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carsul, 'S ', 'S ', carsul), cal(p['0_8'][0], p['C_S_S_C_1'][0], i), cal(p['0_8'][1], p['C_S_S_C_1'][1], i), cal(p['0_8'][2], p['C_S_S_C_1'][2], i), cal(p['0_8'][3], p['C_S_S_C_1'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carsul, inthyd), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carsul, inthyd), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carsul, 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carsul, 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), intsul, cal(p['0_S'][2], p['SH'][2], i), cal(p['0_S'][3], p['SH'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), carsul, cal(p['SH'][2], p['CT'][2], i), cal(p['SH'][3], p['CT'][3], i))

def stock_add_to_all_N(var=variablemake()):
	intsul = var[0]
	newhyd = var[1]
	carsul = var[2]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(intsul, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(newhyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(carsul, 'S', 'sp3')
	p = {}
	with open('Param_files/Stock/Stock.param', 'r') as b:
		data = b.readlines()[1:]
	for line in data:
		p[line.split()[0]] = []
		for point in line.split()[1:]:
			p[line.split()[0]].append(float(point))
	b.close()
	for i in range(11):
		a = i*10
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), intsul, cal(p['0_S'][0], p['SH'][0], i), cal(p['0_S'][1], p['SH'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), carsul, cal(p['SH'][0], p['CT'][0], i), cal(p['SH'][1], p['CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', intsul), cal(p['CT_mS'][0], p['CT_S'][0], i), cal(p['CT_mS'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(intsul, carsul), cal(p['CT_mS'][0], p['CT_S'][0], i), cal(p['CT_mS'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carsul, newhyd), cal(p['CT_mSH'][0], p['CT_HC'][0], i), cal(p['CT_mSH'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carsul, 'S '), cal(p['S_S'][0], p['CT_S'][0], i), cal(p['S_S'][1], p['CT_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', intsul), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', intsul), cal(p['C_C_S'][0], p['C_C_S'][0], i), cal(p['C_C_S'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intsul, carsul), cal(p['Dritt'][0], p['C_S_C'][0], i), cal(p['Dritt'][1], p['C_S_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intsul, carsul, newhyd), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, carsul, newhyd), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intsul, carsul, 'S '), cal(p['C_S_S'][0], p['S_C_S'][0], i), cal(p['C_S_S'][1], p['S_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, carsul, 'S '), cal(p['C_S_S'][0], p['C_C_H'][0], i), cal(p['C_S_S'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carsul, 'S ', 'CT'), cal(p['C_S_S'][0], p['C_S_S'][0], i), cal(p['C_S_S'][1], p['C_S_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carsul, 'S ', '2C'), cal(p['C_S_S'][0], p['C_S_S'][0], i), cal(p['C_S_S'][1], p['C_S_S'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'S ', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'S ', '2C'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'S ', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'S ', '2C'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', intsul, carsul), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', intsul, carsul), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intsul, carsul, 'S '), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intsul, carsul, newhyd), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carsul, 'S ', 'S ', carsul), cal(p['0_1'][0], p['C_S_S_C_2'][0], i), cal(p['0_1'][1], p['C_S_S_C_2'][1], i), cal(p['0_1'][2], p['C_S_S_C_2'][2], i), cal(p['0_1'][3], p['C_S_S_C_2'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carsul, 'S ', 'S ', carsul), cal(p['0_8'][0], p['C_S_S_C_1'][0], i), cal(p['0_8'][1], p['C_S_S_C_1'][1], i), cal(p['0_8'][2], p['C_S_S_C_1'][2], i), cal(p['0_8'][3], p['C_S_S_C_1'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carsul, inthyd), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carsul, inthyd), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carsul, 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carsul, 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), intsul, cal(p['0_S'][2], p['SH'][2], i), cal(p['0_S'][3], p['SH'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), carsul, cal(p['SH'][2], p['CT'][2], i), cal(p['SH'][3], p['CT'][3], i))

def stock_add_to_all_C(var=variablemake()):
	intsul = var[0]
	newhyd = var[1]
	carsul = var[2]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(intsul, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(newhyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(carsul, 'S', 'sp3')
	p = {}
	with open('Param_files/Stock/Stock.param', 'r') as b:
		data = b.readlines()[1:]
	for line in data:
		p[line.split()[0]] = []
		for point in line.split()[1:]:
			p[line.split()[0]].append(float(point))
	b.close()
	for i in range(11):
		a = i*10
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), intsul, cal(p['0_S'][0], p['SH'][0], i), cal(p['0_S'][1], p['SH'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), carsul, cal(p['SH'][0], p['CT'][0], i), cal(p['SH'][1], p['CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', intsul), cal(p['CT_mS'][0], p['CT_S'][0], i), cal(p['CT_mS'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(intsul, carsul), cal(p['CT_mS'][0], p['CT_S'][0], i), cal(p['CT_mS'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carsul, newhyd), cal(p['CT_mSH'][0], p['CT_HC'][0], i), cal(p['CT_mSH'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carsul, 'S '), cal(p['S_S'][0], p['CT_S'][0], i), cal(p['S_S'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carsul, 'CT'), cal(p['CT_S'][0], p['CT_CT'][0], i), cal(p['CT_S'][1], p['CT_CT'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', intsul), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', carsul), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', intsul), cal(p['C_C_S'][0], p['C_C_S'][0], i), cal(p['C_C_S'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', 'CT', carsul), cal(p['S_C_S'][0], p['C_C_S'][0], i), cal(p['S_C_S'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intsul, carsul), cal(p['Dritt'][0], p['C_S_C'][0], i), cal(p['Dritt'][1], p['C_S_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', carsul, intsul), cal(p['C_S_C'][0], p['C_C_S'][0], i), cal(p['C_S_C'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intsul, carsul, newhyd), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, carsul, newhyd), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intsul, carsul, 'S '), cal(p['C_S_S'][0], p['S_C_S'][0], i), cal(p['C_S_S'][1], p['S_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, carsul, 'S '), cal(p['C_S_S'][0], p['C_C_H'][0], i), cal(p['C_S_S'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, carsul, 'CT'), cal(p['C_C_S'][0], p['C_C_H'][0], i), cal(p['C_C_S'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carsul, 'S ', 'CT'), cal(p['C_S_S'][0], p['C_S_S'][0], i), cal(p['C_S_S'][1], p['C_S_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carsul, 'S ', '2C'), cal(p['C_S_S'][0], p['C_S_S'][0], i), cal(p['C_S_S'][1], p['C_S_S'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'S ', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'S ', '2C'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'CT', 'H1'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(intsul, carsul, 'CT', 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'CT', 'H1'), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'CT', 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'S ', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, carsul, 'S ', '2C'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', intsul, carsul), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', intsul, carsul), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intsul, carsul, 'S '), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intsul, carsul, 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intsul, carsul, newhyd), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carsul, 'S ', 'S ', carsul), cal(p['0_1'][0], p['C_S_S_C_2'][0], i), cal(p['0_1'][1], p['C_S_S_C_2'][1], i), cal(p['0_1'][2], p['C_S_S_C_2'][2], i), cal(p['0_1'][3], p['C_S_S_C_2'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carsul, 'S ', 'S ', carsul), cal(p['0_8'][0], p['C_S_S_C_1'][0], i), cal(p['0_8'][1], p['C_S_S_C_1'][1], i), cal(p['0_8'][2], p['C_S_S_C_1'][2], i), cal(p['0_8'][3], p['C_S_S_C_1'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carsul, inthyd), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carsul, inthyd), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carsul, 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
#		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carsul, 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), intsul, cal(p['0_S'][2], p['SH'][2], i), cal(p['0_S'][3], p['SH'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), carsul, cal(p['SH'][2], p['CT'][2], i), cal(p['SH'][3], p['CT'][3], i))
