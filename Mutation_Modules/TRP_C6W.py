# TRP to C6W Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/TRP.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/C6W.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                #changeLJPair(parm, ':{}@HH2 :{}@HZ1 0 0'.format(vxi, vxi)).execute()
                #changeLJPair(parm, ':{}@HE3 :{}@HZ2 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD1'.format(vxi), bc['CD1']+((fc['CD1']-bc['CD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD1'.format(vxi), bc['HD1']+((fc['HD1']-bc['HD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@NE1'.format(vxi), bc['NE1']+((fc['NE1']-bc['NE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), bc['HE1']+((fc['HE1']-bc['HE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE2'.format(vxi), bc['CE2']+((fc['CE2']-bc['CE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ2'.format(vxi), bc['CZ2']+((fc['CZ2']-bc['CZ2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vxi), bc['HZ2']+((fc['HZ2']-bc['HZ2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CH2'.format(vxi), bc['CH2']+((fc['CH2']-bc['CH2'])/10)*i).execute()
                change(parm, 'charge', ':{}@ClH2'.format(vxi), bc['HH2']+((fc['ClH2']-bc['HH2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ3'.format(vxi), bc['CZ3']+((fc['CZ3']-bc['CZ3'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ3'.format(vxi), bc['HZ3']+((fc['HZ3']-bc['HZ3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE3'.format(vxi), bc['CE3']+((fc['CE3']-bc['CE3'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), bc['HE3']+((fc['HE3']-bc['HE3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']+((fc['CD2']-bc['CD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'HH2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('ClH2'))
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

def makevxi_alt(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	CD1 = struct.residue_dict[aa].atom_dict['CD1']
	CD2 = struct.residue_dict[aa].atom_dict['CD2']
	HD1 = struct.residue_dict[aa].atom_dict['HD1']
	CE = struct.residue_dict[aa].atom_dict['CE']
	HE = struct.residue_dict[aa].atom_dict['HE']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'CD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CD2'))
			elif atom.get_name() == 'HD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CE3'))
			elif atom.get_name() == 'CD2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CD1'))
				pdb.write(atom.halfway_between('NE1', CD2, CE))
				pdb.write(atom.superimposed1('HE1', CD2))
			elif atom.get_name() == 'HD2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HD1'))
				pdb.write(atom.superimposed1('HE3', CD1))
			elif atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CE2'))
			elif atom.get_name() == 'HE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ2'))
				pdb.write(atom.superimposed1('HZ2', CE))
				pdb.write(atom.thirdone_between('CZ3', HD1, HE))
				pdb.write(atom.superimposed1('HZ3', HD1))
				pdb.write(atom.thirdtwo_between('CH2', HD1, HE))
				pdb.write(atom.superimposed1('HH2', HE))
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
	var16 = sym + 'f'
	var17 = sym + 'g'
	var18 = sym + 'h'
	return var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, var15, var16, var17, var18

def lib_make(ff, outputfile, vxi='VXI', var=variablemake()):
	chloro = var[0]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/TRP-C6W.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "N"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "Cl"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "C"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "C"\n'%vxi)
	ctrl.write('set %s.1.23 element "C"\n'%vxi)
	ctrl.write('set %s.1.24 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "CD1"\n'%vxi)
	ctrl.write('set %s.1.10 name "HD1"\n'%vxi)
	ctrl.write('set %s.1.11 name "NE1"\n'%vxi)
	ctrl.write('set %s.1.12 name "HE1"\n'%vxi)
	ctrl.write('set %s.1.13 name "CE2"\n'%vxi)
	ctrl.write('set %s.1.14 name "CZ2"\n'%vxi)
	ctrl.write('set %s.1.15 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.16 name "CH2"\n'%vxi)
	ctrl.write('set %s.1.17 name "ClH2"\n'%vxi)
	ctrl.write('set %s.1.18 name "CZ3"\n'%vxi)
	ctrl.write('set %s.1.19 name "HZ3"\n'%vxi)
	ctrl.write('set %s.1.20 name "CE3"\n'%vxi)
	ctrl.write('set %s.1.21 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.22 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.23 name "C"\n'%vxi)
	ctrl.write('set %s.1.24 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "C*"\n'%vxi)
	ctrl.write('set %s.1.9 type "CW"\n'%vxi)
	ctrl.write('set %s.1.10 type "H4"\n'%vxi)
	ctrl.write('set %s.1.11 type "NA"\n'%vxi)
	ctrl.write('set %s.1.12 type "H"\n'%vxi)
	ctrl.write('set %s.1.13 type "CN"\n'%vxi)
	ctrl.write('set %s.1.14 type "CA"\n'%vxi)
	ctrl.write('set %s.1.15 type "HA"\n'%vxi)
	ctrl.write('set %s.1.16 type "CA"\n'%vxi)
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, chloro))
	ctrl.write('set %s.1.18 type "CA"\n'%vxi)
	ctrl.write('set %s.1.19 type "HA"\n'%vxi)
	ctrl.write('set %s.1.20 type "CA"\n'%vxi)
	ctrl.write('set %s.1.21 type "HA"\n'%vxi)
	ctrl.write('set %s.1.22 type "CB"\n'%vxi)
	ctrl.write('set %s.1.23 type "C"\n'%vxi)
	ctrl.write('set %s.1.24 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
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

def lac(x, y, i):
	num = y+((x-y)/10)*i
	return num

def stock_add_to_all(var=variablemake()):
	chloro = var[0]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(chloro, 'Cl', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), chloro, cal(p['HA'][0], p['Cl'][0], i), cal(p['HA'][1], p['Cl'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CA', chloro), cal(p['CA_HA'][0], p['CA_Cl'][0], i), cal(p['CA_HA'][1], p['CA_Cl'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', 'CA', chloro), cal(p['F_CA_CA_HA'][0], p['CA_CA_Cl'][0], i), cal(p['F_CA_CA_HA'][1], p['CA_CA_Cl'][1], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', 'CA', chloro), cal(p['Imp_0'][0], p['Ring_imp'][0], i), cal(p['Imp_0'][1], p['Ring_imp'][1], i), cal(p['Imp_0'][2], p['Ring_imp'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), chloro, cal(p['HA'][2], p['Cl'][2], i), cal(p['HA'][3], p['Cl'][3], i))
