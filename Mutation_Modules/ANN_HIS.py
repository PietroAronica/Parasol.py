# ANN to HIS Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
	bc = {}
        with open('Param_files/AminoAcid/ANN.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/HIS.param', 'r') as b:
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
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@ND1'.format(vxi), (fc['ND1']/10)*i).execute()
                change(parm, 'charge', ':{}@CE1'.format(vxi), bc['CD1']+((fc['CE1']-bc['CD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), bc['HD1']+((fc['HE1']-bc['HD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@NE2'.format(vxi), bc['CE']+((fc['NE2']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), bc['HE']+((fc['HE2']-bc['HE'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']+((fc['CD2']-bc['CD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	CG = struct.residue_dict[aa].atom_dict['CG']
	CD1 = struct.residue_dict[aa].atom_dict['CD1']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'CG' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('ND1', CG, CD1)) 
			elif atom.get_name() == 'CD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CE1')) 
			elif atom.get_name() == 'HD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HE1')) 
			elif atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('NE2')) 
			elif atom.get_name() == 'HE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HE2')) 
			else:
				pdb.write(atom.formatted())
                try:
                        pdb.write(struct.other_dict[res.get_resnumber()].ter())
                except:
                        pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def lib_make(ff, outputfile, vxi='VXI', cg='1c', nd1='2n', ce1='3c', he1='3h', ne2='4n', he2='4h', cd2='5c', hd2='5h'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ANN-HIS.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "N"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "N"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "ND1"\n'%vxi)
	ctrl.write('set %s.1.10 name "CE1"\n'%vxi)
	ctrl.write('set %s.1.11 name "HE1"\n'%vxi)
	ctrl.write('set %s.1.12 name "NE2"\n'%vxi)
	ctrl.write('set %s.1.13 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.14 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.15 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.16 name "C"\n'%vxi)
	ctrl.write('set %s.1.17 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, cg))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, nd1))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, ce1))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, he1))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, ne2))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, he2))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, cd2))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, hd2))
	ctrl.write('set %s.1.16 type "C"\n'%vxi)
	ctrl.write('set %s.1.17 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
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

def stock_add_to_all(cg='1c', nd1='2n', ce1='3c', he1='3h', ne2='4n', he2='4h', cd2='5c', hd2='5h'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(cg, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(nd1, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(ce1, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(he1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ne2, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(he2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cd2, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hd2, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cg, cal(p['CA'][0], p['CA'][0], i), cal(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), nd1, cal(p['0_N'][0], p['NA'][0], i), cal(p['0_N'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ce1, cal(p['CA'][0], p['CA'][0], i), cal(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he1, cal(p['HA'][0], p['H5'][0], i), cal(p['HA'][1], p['H5'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ne2, cal(p['CA'][0], p['NA'][0], i), cal(p['CA'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he2, cal(p['HA'][0], p['H'][0], i), cal(p['HA'][1], p['H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd2, cal(p['CA'][0], p['CA'][0], i), cal(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd2, cal(p['HA'][0], p['H4'][0], i), cal(p['HA'][1], p['H4'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', cg), cal(p['CT_CA'][0], p['CT_CC'][0], i), cal(p['CT_CA'][1], p['CT_CC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cg, nd1), cal(p['NB_mCA1'][0], p['CC_NB'][0], i), cal(p['NB_mCA1'][1], p['CC_NB'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(nd1, ce1), cal(p['NB_mCA2'][0], p['CR_NB'][0], i), cal(p['NB_mCA2'][1], p['CR_NB'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce1, he1), cal(p['CA_HA'][0], p['CA_HA'][0], i), cal(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce1, ne2), cal(p['CA_CA'][0], p['CR_NA'][0], i), cal(p['CA_CA'][1], p['CR_NA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne2, he2), cal(p['CA_HA'][0], p['NA_H'][0], i), cal(p['CA_HA'][1], p['NA_H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne2, cd2), cal(p['CA_CA'][0], p['CW_NA'][0], i), cal(p['CA_CA'][1], p['CW_NA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd2, hd2), cal(p['CA_HA'][0], p['CA_HA'][0], i), cal(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd2, cg), cal(p['CA_CA'][0], p['CC_CW'][0], i), cal(p['CA_CA'][1], p['CC_CW'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', cg), cal(p['C_C_CA'][0], p['C_C_CC'][0], i), cal(p['C_C_CA'][1], p['C_C_CC'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', cg), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', cg, nd1), cal(p['A_CT_CA_CA'][0], p['F_CT_CA_CA'][0], i), cal(p['A_CT_CA_CA'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', cg, cd2), cal(p['A_CT_CA_CA'][0], p['F_CT_CA_CA'][0], i), cal(p['A_CT_CA_CA'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(nd1, cg, cd2), cal(p['A_CA_CA_CA'][0], p['F_CT_CA_CA'][0], i), cal(p['A_CA_CA_CA'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cg, nd1, ce1), cal(p['Dritt'][0], p['C_C_O2'][0], i), cal(p['Dritt'][1], p['C_C_O2'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(nd1, ce1, ne2), cal(p['A_CA_CA_CA'][0], p['F_CT_CA_CA'][0], i), cal(p['A_CA_CA_CA'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(nd1, ce1, he1), cal(p['A_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['A_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne2, ce1, he1), cal(p['A_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['A_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce1, ne2, cd2), cal(p['A_CA_CA_CA'][0], p['F_CT_CA_CA'][0], i), cal(p['A_CA_CA_CA'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce1, ne2, he2), cal(p['A_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['A_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd2, ne2, he2), cal(p['A_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['A_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne2, cd2, cg), cal(p['A_CA_CA_CA'][0], p['F_CT_CA_CA'][0], i), cal(p['A_CA_CA_CA'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne2, cd2, hd2), cal(p['A_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['A_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cg, cd2, hd2), cal(p['A_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['A_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, nd1), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, cd2), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, nd1), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd2), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, nd1, ce1), cal(p['0_6'][0], p['X_CC_NB_X'][0], i), cal(p['0_6'][1], p['X_CC_NB_X'][1], i), cal(p['0_6'][2], p['X_CC_NB_X'][2], i), cal(p['0_6'][3], p['X_CC_NB_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, cg, nd1, ce1), cal(p['0_6'][0], p['X_CC_NB_X'][0], i), cal(p['0_6'][1], p['X_CC_NB_X'][1], i), cal(p['0_6'][2], p['X_CC_NB_X'][2], i), cal(p['0_6'][3], p['X_CC_NB_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, ne2), cal(p['Ring_0'][0], p['X_CC_CW_X'][0], i), cal(p['Ring_0'][1], p['X_CC_CW_X'][1], i), cal(p['Ring_0'][2], p['X_CC_CW_X'][2], i), cal(p['Ring_0'][3], p['X_CC_CW_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, hd2), cal(p['Ring_0'][0], p['X_CC_CW_X'][0], i), cal(p['Ring_0'][1], p['X_CC_CW_X'][1], i), cal(p['Ring_0'][2], p['X_CC_CW_X'][2], i), cal(p['Ring_0'][3], p['X_CC_CW_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nd1, cg, cd2, ne2), cal(p['Ring_0'][0], p['X_CC_CW_X'][0], i), cal(p['Ring_0'][1], p['X_CC_CW_X'][1], i), cal(p['Ring_0'][2], p['X_CC_CW_X'][2], i), cal(p['Ring_0'][3], p['X_CC_CW_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nd1, cg, cd2, hd2), cal(p['Ring_0'][0], p['X_CC_CW_X'][0], i), cal(p['Ring_0'][1], p['X_CC_CW_X'][1], i), cal(p['Ring_0'][2], p['X_CC_CW_X'][2], i), cal(p['Ring_0'][3], p['X_CC_CW_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, ne2, he2), cal(p['Ring_0'][0], p['X_CW_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CW_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CW_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CW_NA_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, ne2, ce1), cal(p['Ring_0'][0], p['X_CW_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CW_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CW_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CW_NA_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd2, cd2, ne2, he2), cal(p['Ring_0'][0], p['X_CW_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CW_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CW_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CW_NA_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd2, cd2, ne2, ce1), cal(p['Ring_0'][0], p['X_CW_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CW_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CW_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CW_NA_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, nd1, ce1, he1), cal(p['0_6'][0], p['X_CR_NB_X'][0], i), cal(p['0_6'][1], p['X_CR_NB_X'][1], i), cal(p['0_6'][2], p['X_CR_NB_X'][2], i), cal(p['0_6'][3], p['X_CR_NB_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, nd1, ce1, ne2), cal(p['0_6'][0], p['X_CR_NB_X'][0], i), cal(p['0_6'][1], p['X_CR_NB_X'][1], i), cal(p['0_6'][2], p['X_CR_NB_X'][2], i), cal(p['0_6'][3], p['X_CR_NB_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nd1, ce1, ne2, he2), cal(p['Ring_0'][0], p['X_CR_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CR_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CR_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CR_NA_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nd1, ce1, ne2, cd2), cal(p['Ring_0'][0], p['X_CR_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CR_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CR_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CR_NA_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he1, ce1, ne2, he2), cal(p['Ring_0'][0], p['X_CR_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CR_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CR_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CR_NA_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he1, ce1, ne2, cd2), cal(p['Ring_0'][0], p['X_CR_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CR_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CR_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CR_NA_X'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne2, cd2, cg, hd2), cal(p['Ring_imp'][0], p['Ring_imp'][0], i), cal(p['Ring_imp'][1], p['Ring_imp'][1], i), cal(p['Ring_imp'][2], p['Ring_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne2, ce1, nd1, he1), cal(p['Ring_imp'][0], p['Ring_imp'][0], i), cal(p['Ring_imp'][1], p['Ring_imp'][1], i), cal(p['Ring_imp'][2], p['Ring_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cd2, cg, nd1), cal(p['Ring_imp'][0], p['Ring_imp'][0], i), cal(p['Ring_imp'][1], p['Ring_imp'][1], i), cal(p['Ring_imp'][2], p['Ring_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cd2, cg, ne2), cal(p['Ring_imp'][0], p['Ring_imp'][0], i), cal(p['Ring_imp'][1], p['Ring_imp'][1], i), cal(p['Ring_imp'][2], p['Ring_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, ne2, ce1, he2), cal(p['Ring_imp'][0], p['Ami_imp'][0], i), cal(p['Ring_imp'][1], p['Ami_imp'][1], i), cal(p['Ring_imp'][2], p['Ami_imp'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cg, cal(p['CA'][2], p['CA'][2], i), cal(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), nd1, cal(p['0_N'][2], p['NA'][2], i), cal(p['0_N'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ce1, cal(p['CA'][2], p['CA'][2], i), cal(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he1, cal(p['HA'][2], p['H5'][2], i), cal(p['HA'][3], p['H5'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ne2, cal(p['CA'][2], p['NA'][2], i), cal(p['CA'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he2, cal(p['HA'][2], p['H'][2], i), cal(p['HA'][3], p['H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd2, cal(p['CA'][2], p['CA'][2], i), cal(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd2, cal(p['HA'][2], p['H4'][2], i), cal(p['HA'][3], p['H4'][3], i))
