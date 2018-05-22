# LEU to B1L Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/LEU.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/B1L.param', 'r') as b:
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
                change(parm, 'charge', ':{}@CA1'.format(vxi), bc['CA']+((fc['CA1']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA1'.format(vxi), bc['HA']+((fc['HA1']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG'.format(vxi), bc['HG']+((fc['HG']-bc['HG'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD1'.format(vxi), bc['CD1']+((fc['CD1']-bc['CD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD11'.format(vxi), bc['HD11']+((fc['HD11']-bc['HD11'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD12'.format(vxi), bc['HD12']+((fc['HD12']-bc['HD12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD13'.format(vxi), bc['HD13']+((fc['HD13']-bc['HD13'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']+((fc['CD2']-bc['CD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD21'.format(vxi), bc['HD21']+((fc['HD21']-bc['HD21'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD22'.format(vxi), bc['HD22']+((fc['HD22']-bc['HD22'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD23'.format(vxi), bc['HD23']+((fc['HD23']-bc['HD23'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA2'.format(vxi), (fc['CA2']/10)*i).execute()
                change(parm, 'charge', ':{}@HA22'.format(vxi), (fc['HA22']/10)*i).execute()
                change(parm, 'charge', ':{}@HA23'.format(vxi), (fc['HA23']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CA = struct.residue_dict[aa].atom_dict['CA']
        C = struct.residue_dict[aa].atom_dict['C']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'CA' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CA1')) 
			elif atom.get_name() == 'HA' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HA1')) 
			elif atom.get_name() == 'HD23' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('CA2', CA, C)) 
				pdb.write(atom.superimposed1('HA22', C)) 
				pdb.write(atom.superimposed2('HA23', C)) 
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

def lib_make(ff, outputfile, vxi='VXI', ic='1c', ih='1h'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/LEU-B1L.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
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
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "C"\n'%vxi)
	ctrl.write('set %s.1.22 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA1"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA1"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG"\n'%vxi)
	ctrl.write('set %s.1.10 name "CD1"\n'%vxi)
	ctrl.write('set %s.1.11 name "HD11"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD12"\n'%vxi)
	ctrl.write('set %s.1.13 name "HD13"\n'%vxi)
	ctrl.write('set %s.1.14 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.15 name "HD21"\n'%vxi)
	ctrl.write('set %s.1.16 name "HD22"\n'%vxi)
	ctrl.write('set %s.1.17 name "HD23"\n'%vxi)
	ctrl.write('set %s.1.18 name "CA2"\n'%vxi)
	ctrl.write('set %s.1.19 name "HA22"\n'%vxi)
	ctrl.write('set %s.1.20 name "HA23"\n'%vxi)
	ctrl.write('set %s.1.21 name "C"\n'%vxi)
	ctrl.write('set %s.1.22 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
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
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, ic))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, ih))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, ih))
	ctrl.write('set %s.1.21 type "C"\n'%vxi)
	ctrl.write('set %s.1.22 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
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
	ctrl.write('bond %s.1.18 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.21\n'%(vxi, vxi))
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

def stock_add_to_all(ic='1c', ih='1h'):
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
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ic, 'C '), cal(p['Halve'][0], p['CT_C'][0], i), cal(p['Halve'][1], p['CT_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N ', 'CT', ic), cal(p['N_CT_C'][0], p['CT_CT_N'][0], i), cal(p['N_CT_C'][1], p['CT_CT_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', ic), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', ic), cal(p['CT_CT_C'][0], p['C_C_C'][0], i), cal(p['CT_CT_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ic, ih), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ic, 'C '), cal(p['Dritt'][0], p['CT_CT_C'][0], i), cal(p['Dritt'][1], p['CT_CT_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ih, ic, 'C '), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ih, ic, ih), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ic, 'C ', 'O '), cal(p['C_C_O'][0], p['C_C_O'][0], i), cal(p['C_C_O'][1], p['C_C_O'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ic, 'C ', 'N '), cal(p['C_C_N'][0], p['C_C_N'][0], i), cal(p['C_C_N'][1], p['C_C_N'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', ic, ih), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', ic, 'C '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'CT', ic, ih), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'CT', ic, 'C '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ic, ih), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ic, 'C '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'C ', 'N '), cal(p['0_3'][0], p['C_C_C_N_1'][0], i), cal(p['0_3'][1], p['C_C_C_N_1'][1], i), cal(p['0_3'][2], p['C_C_C_N_1'][2], i), cal(p['0_3'][3], p['C_C_C_N_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'C ', 'N '), cal(p['0_8'][0], p['C_C_C_N_2'][0], i), cal(p['0_8'][1], p['C_C_C_N_2'][1], i), cal(p['0_8'][2], p['C_C_C_N_2'][2], i), cal(p['0_8'][3], p['C_C_C_N_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'C ', 'N '), cal(p['0_2'][0], p['C_C_C_N_3'][0], i), cal(p['0_2'][1], p['C_C_C_N_3'][1], i), cal(p['0_2'][2], p['C_C_C_N_3'][2], i), cal(p['0_2'][3], p['C_C_C_N_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'C ', 'N '), cal(p['0_7'][0], p['C_C_C_N_4'][0], i), cal(p['0_7'][1], p['C_C_C_N_4'][1], i), cal(p['0_7'][2], p['C_C_C_N_4'][2], i), cal(p['0_7'][3], p['C_C_C_N_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ih, ic, 'C ', 'N '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ic, 'C ', 'O '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ih, ic, 'C ', 'O '), cal(p['0_10'][0], p['H_C_C_O_1'][0], i), cal(p['0_10'][1], p['H_C_C_O_1'][1], i), cal(p['0_10'][2], p['H_C_C_O_1'][2], i), cal(p['0_10'][3], p['H_C_C_O_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ih, ic, 'C ', 'O '), cal(p['0_8'][0], p['H_C_C_O_2'][0], i), cal(p['0_8'][1], p['H_C_C_O_2'][1], i), cal(p['0_8'][2], p['H_C_C_O_2'][2], i), cal(p['0_8'][3], p['H_C_C_O_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ih, ic, 'C ', 'O '), cal(p['0_9'][0], p['H_C_C_O_3'][0], i), cal(p['0_9'][1], p['H_C_C_O_3'][1], i), cal(p['0_9'][2], p['H_C_C_O_3'][2], i), cal(p['0_9'][3], p['H_C_C_O_3'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ic, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ih, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
