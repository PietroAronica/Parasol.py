# VAL to B2V Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/VAL.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/B2V.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA1'.format(vxi), (fc['CA1']/10)*i).execute()
                change(parm, 'charge', ':{}@HA12'.format(vxi), (fc['HA12']/10)*i).execute()
                change(parm, 'charge', ':{}@HA13'.format(vxi), (fc['HA13']/10)*i).execute()
                change(parm, 'charge', ':{}@CA2'.format(vxi), bc['CA']+((fc['CA2']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA2'.format(vxi), bc['HA']+((fc['HA2']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB'.format(vxi), bc['HB']+((fc['HB']-bc['HB'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG1'.format(vxi), bc['CG1']+((fc['CG1']-bc['CG1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG11'.format(vxi), bc['HG11']+((fc['HG11']-bc['HG11'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG12'.format(vxi), bc['HG12']+((fc['HG12']-bc['HG12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG12'.format(vxi), bc['HG13']+((fc['HG13']-bc['HG13'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG2'.format(vxi), bc['CG2']+((fc['CG2']-bc['CG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG21'.format(vxi), bc['HG21']+((fc['HG21']-bc['HG21'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG22'.format(vxi), bc['HG22']+((fc['HG22']-bc['HG22'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG22'.format(vxi), bc['HG23']+((fc['HG23']-bc['HG23'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CA = struct.residue_dict[aa].atom_dict['CA']
        N = struct.residue_dict[aa].atom_dict['N']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'H' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('CA1', CA, N)) 
				pdb.write(atom.superimposed1('HA12', N)) 
				pdb.write(atom.superimposed2('HA13', N)) 
			elif atom.get_name() == 'CA' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CA2')) 
			elif atom.get_name() == 'HA' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HA2')) 
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
	ic = var[0]
	ih = var[1]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/VAL-B2V.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "H"\n'%vxi)
	ctrl.write('set %s.1.6 element "C"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA1"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA12"\n'%vxi)
	ctrl.write('set %s.1.5 name "HA13"\n'%vxi)
	ctrl.write('set %s.1.6 name "CA2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HA2"\n'%vxi)
	ctrl.write('set %s.1.8 name "CB"\n'%vxi)
	ctrl.write('set %s.1.9 name "HB"\n'%vxi)
	ctrl.write('set %s.1.10 name "CG1"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG11"\n'%vxi)
	ctrl.write('set %s.1.12 name "HG12"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG13"\n'%vxi)
	ctrl.write('set %s.1.14 name "CG2"\n'%vxi)
	ctrl.write('set %s.1.15 name "HG21"\n'%vxi)
	ctrl.write('set %s.1.16 name "HG22"\n'%vxi)
	ctrl.write('set %s.1.17 name "HG23"\n'%vxi)
	ctrl.write('set %s.1.18 name "C"\n'%vxi)
	ctrl.write('set %s.1.19 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "%s"\n'%(vxi, ic))
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, ih))
	ctrl.write('set %s.1.5 type "%s"\n'%(vxi, ih))
	ctrl.write('set %s.1.6 type "CT"\n'%vxi)
	ctrl.write('set %s.1.7 type "H1"\n'%vxi)
	ctrl.write('set %s.1.8 type "CT"\n'%vxi)
	ctrl.write('set %s.1.9 type "HC"\n'%vxi)
	ctrl.write('set %s.1.10 type "CT"\n'%vxi)
	ctrl.write('set %s.1.11 type "HC"\n'%vxi)
	ctrl.write('set %s.1.12 type "HC"\n'%vxi)
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "CT"\n'%vxi)
	ctrl.write('set %s.1.15 type "HC"\n'%vxi)
	ctrl.write('set %s.1.16 type "HC"\n'%vxi)
	ctrl.write('set %s.1.17 type "HC"\n'%vxi)
	ctrl.write('set %s.1.18 type "C"\n'%vxi)
	ctrl.write('set %s.1.19 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
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
	ic = var[0]
	ih = var[1]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(ic, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(ih, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ic, cal(p['0_C'][0], p['CT'][0], i), cal(p['0_C'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ih, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', ic), cal(p['Halve'][0], p['CT_CT'][0], i), cal(p['Halve'][1], p['CT_CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ic, ih), cal(p['Halve'][0], p['CT_HC'][0], i), cal(p['Halve'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ic, 'N '), cal(p['Halve'][0], p['CT_N'][0], i), cal(p['Halve'][1], p['CT_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', 'CT', ic), cal(p['N_CT_C'][0], p['CT_CT_C'][0], i), cal(p['N_CT_C'][1], p['CT_CT_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', ic), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', ic), cal(p['CT_CT_N'][0], p['C_C_C'][0], i), cal(p['CT_CT_N'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ic, ih), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ic, 'N '), cal(p['Dritt'][0], p['CT_CT_N'][0], i), cal(p['Dritt'][1], p['CT_CT_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ih, ic, 'N '), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ih, ic, ih), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ic, 'N ', 'H '), cal(p['CT_N_H'][0], p['CT_N_H'][0], i), cal(p['CT_N_H'][1], p['CT_N_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ic, 'N ', 'C '), cal(p['C_N_CT'][0], p['C_N_CT'][0], i), cal(p['C_N_CT'][1], p['C_N_CT'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', ic, ih), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', ic, 'N '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'CT', ic, ih), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'CT', ic, 'N '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ic, ih), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ic, 'N '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'N ', 'C '), cal(p['0_3'][0], p['C_C_N_C_1'][0], i), cal(p['0_3'][1], p['C_C_N_C_1'][1], i), cal(p['0_3'][2], p['C_C_N_C_1'][2], i), cal(p['0_3'][3], p['C_C_N_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'N ', 'C '), cal(p['0_8'][0], p['C_C_N_C_2'][0], i), cal(p['0_8'][1], p['C_C_N_C_2'][1], i), cal(p['0_8'][2], p['C_C_N_C_2'][2], i), cal(p['0_8'][3], p['C_C_N_C_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'N ', 'C '), cal(p['0_2'][0], p['C_C_N_C_3'][0], i), cal(p['0_2'][1], p['C_C_N_C_3'][1], i), cal(p['0_2'][2], p['C_C_N_C_3'][2], i), cal(p['0_2'][3], p['C_C_N_C_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'N ', 'C '), cal(p['0_7'][0], p['C_C_N_C_4'][0], i), cal(p['0_7'][1], p['C_C_N_C_4'][1], i), cal(p['0_7'][2], p['C_C_N_C_4'][2], i), cal(p['0_7'][3], p['C_C_N_C_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ih, ic, 'N ', 'C '), cal(p['X_C_N_X'][0], p['X_C_N_X'][0], i), cal(p['X_C_N_X'][1], p['X_C_N_X'][1], i), cal(p['X_C_N_X'][2], p['X_C_N_X'][2], i), cal(p['X_C_N_X'][3], p['X_C_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'N ', 'H '), cal(p['X_C_N_X'][0], p['X_C_N_X'][0], i), cal(p['X_C_N_X'][1], p['X_C_N_X'][1], i), cal(p['X_C_N_X'][2], p['X_C_N_X'][2], i), cal(p['X_C_N_X'][3], p['X_C_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ih, ic, 'N ', 'H '), cal(p['X_C_N_X'][0], p['X_C_N_X'][0], i), cal(p['X_C_N_X'][1], p['X_C_N_X'][1], i), cal(p['X_C_N_X'][2], p['X_C_N_X'][2], i), cal(p['X_C_N_X'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ic, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ih, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
