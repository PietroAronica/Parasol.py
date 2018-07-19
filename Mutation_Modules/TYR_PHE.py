# TYR to PHE Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/TYR.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/PHE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HZ :{}@HH 0 0'.format(vxi, vxi)).execute()
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
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']+((fc['CD2']-bc['CD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE1'.format(vxi), bc['CE1']+((fc['CE1']-bc['CE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), bc['HE1']+((fc['HE1']-bc['HE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE2'.format(vxi), bc['CE2']+((fc['CE2']-bc['CE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), bc['HE2']+((fc['HE2']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ'.format(vxi), bc['CZ']+((fc['CZ']-bc['CZ'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ'.format(vxi), (fc['HZ']/10)*i).execute()
                change(parm, 'charge', ':{}@OH'.format(vxi), bc['OH']-(bc['OH']/10)*i).execute()
                change(parm, 'charge', ':{}@HH'.format(vxi), bc['HH']-(bc['HH']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	CZ = struct.residue_dict[aa].atom_dict['CZ']
	OH = struct.residue_dict[aa].atom_dict['OH']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'CZ' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.superimposed1('HZ', OH)) 
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
	hz = var[0]
	oh = var[1]
	hh = var[2]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/PHE-TYR.pdb\n"%vxi)
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
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "O"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "C"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "C"\n'%vxi)
	ctrl.write('set %s.1.22 element "O"\n'%vxi)
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
	ctrl.write('set %s.1.11 name "CE1"\n'%vxi)
	ctrl.write('set %s.1.12 name "HE1"\n'%vxi)
	ctrl.write('set %s.1.13 name "CZ"\n'%vxi)
	ctrl.write('set %s.1.14 name "HZ"\n'%vxi)
	ctrl.write('set %s.1.15 name "OH"\n'%vxi)
	ctrl.write('set %s.1.16 name "HH"\n'%vxi)
	ctrl.write('set %s.1.17 name "CE2"\n'%vxi)
	ctrl.write('set %s.1.18 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.19 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.20 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.21 name "C"\n'%vxi)
	ctrl.write('set %s.1.22 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "CA"\n'%vxi)
	ctrl.write('set %s.1.9 type "CA"\n'%vxi)
	ctrl.write('set %s.1.10 type "HA"\n'%vxi)
	ctrl.write('set %s.1.11 type "CA"\n'%vxi)
	ctrl.write('set %s.1.12 type "HA"\n'%vxi)
	ctrl.write('set %s.1.13 type "CA"\n'%vxi)
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, hz))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, oh))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, hh))
	ctrl.write('set %s.1.17 type "CA"\n'%vxi)
	ctrl.write('set %s.1.18 type "HA"\n'%vxi)
	ctrl.write('set %s.1.19 type "CA"\n'%vxi)
	ctrl.write('set %s.1.20 type "HA"\n'%vxi)
	ctrl.write('set %s.1.21 type "C"\n'%vxi)
	ctrl.write('set %s.1.22 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.22\n'%(vxi, vxi))
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
	hz = var[0]
	oh = var[1]
	hh = var[2]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(hz, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(oh, 'O', 'sp3')
	Frcmod_creator.TYPE_insert(hh, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz, lac(p['HA'][0], p['0_H'][0], i), lac(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), oh, lac(p['0_O'][0], p['OH'][0], i), lac(p['0_O'][1], p['OH'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh, lac(p['0_H'][0], p['HO'][0], i), lac(p['0_O'][1], p['HO'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CA', hz), lac(p['CA_HA'][0], p['HA_sOH'][0], i), lac(p['CA_HA'][1], p['HA_sOH'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CA', oh), lac(p['OH_mHA'][0], p['CA_OH'][0], i), lac(p['OH_mHA'][1], p['CA_OH'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(oh, hh), lac(p['HO_mHA'][0], p['OH_HO'][0], i), lac(p['HO_mHA'][1], p['OH_HO'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', 'CA', hz), lac(p['F_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', 'CA', oh), lac(p['F_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz, 'CA', oh), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', oh, hh), lac(p['Dritt'][0], p['CA_O_H'][0], i), lac(p['Dritt'][1], p['CA_O_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz, 'CA', oh, hh), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', oh, hh), lac(p['0_6'][0], p['Ring_Tyr'][0], i), lac(p['0_6'][1], p['Ring_Tyr'][1], i), lac(p['0_6'][2], p['Ring_Tyr'][2], i), lac(p['0_6'][3], p['Ring_Tyr'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', 'CA', oh), lac(p['Ring_0'][0], p['Ring_Dihe'][0], i), lac(p['Ring_0'][1], p['Ring_Dihe'][1], i), lac(p['Ring_0'][2], p['Ring_Dihe'][2], i), lac(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HA', 'CA', 'CA', oh), lac(p['Ring_0'][0], p['Ring_Dihe'][0], i), lac(p['Ring_0'][1], p['Ring_Dihe'][1], i), lac(p['Ring_0'][2], p['Ring_Dihe'][2], i), lac(p['Ring_0'][3], p['Ring_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', 'CA', hz), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HA', 'CA', 'CA', hz), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', 'CA', hz), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', 'CA', oh), lac(p['Imp_0'][0], p['Ring_imp'][0], i), lac(p['Imp_0'][1], p['Ring_imp'][1], i), lac(p['Imp_0'][2], p['Ring_imp'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz, lac(p['HA'][2], p['0_H'][2], i), lac(p['HA'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), oh, lac(p['0_O'][2], p['OH'][2], i), lac(p['0_O'][3], p['OH'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh, lac(p['0_H'][2], p['HO'][2], i), lac(p['0_O'][3], p['HO'][3], i))
