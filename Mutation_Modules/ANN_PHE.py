# ANN to PHE Mutation

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
        with open('Param_files/AminoAcid/PHE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HE1 :{}@HE2 0 0'.format(vxi, vxi)).execute()
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
                change(parm, 'charge', ':{}@CE1'.format(vxi), (fc['CE1']/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), (fc['HE1']/10)*i).execute()
                change(parm, 'charge', ':{}@CE2'.format(vxi), (fc['CE2']/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), (fc['HE2']/10)*i).execute()
                change(parm, 'charge', ':{}@CZ'.format(vxi), bc['CE']+((fc['CZ']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ'.format(vxi), bc['HE']+((fc['HZ']-bc['HE'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	CD1 = struct.residue_dict[aa].atom_dict['CD1']
	CD2 = struct.residue_dict[aa].atom_dict['CD2']
	CE = struct.residue_dict[aa].atom_dict['CE']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'HD2' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('CE1', CD1, CE)) 
				pdb.write(atom.superimposed1('HE1', CE)) 
				pdb.write(atom.halfway_between('CE2', CD2, CE)) 
				pdb.write(atom.superimposed2('HE2', CE)) 
			elif atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ')) 
			elif atom.get_name() == 'HE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HZ')) 
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

def lib_make(ff, outputfile, vxi='VXI', cd='1c', ce='2c', he='2h'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/PHE-ANN.pdb\n"%vxi)
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
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "C"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "O"\n'%vxi)
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
	ctrl.write('set %s.1.11 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.13 name "CE1"\n'%vxi)
	ctrl.write('set %s.1.14 name "HE1"\n'%vxi)
	ctrl.write('set %s.1.15 name "CE2"\n'%vxi)
	ctrl.write('set %s.1.16 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.17 name "CZ"\n'%vxi)
	ctrl.write('set %s.1.18 name "HZ"\n'%vxi)
	ctrl.write('set %s.1.19 name "C"\n'%vxi)
	ctrl.write('set %s.1.20 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "CA"\n'%vxi)
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, cd))
	ctrl.write('set %s.1.10 type "HA"\n'%vxi)
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, cd))
	ctrl.write('set %s.1.12 type "HA"\n'%vxi)
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, ce))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, he))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, ce))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, he))
	ctrl.write('set %s.1.17 type "CA"\n'%vxi)
	ctrl.write('set %s.1.18 type "HA"\n'%vxi)
	ctrl.write('set %s.1.19 type "C"\n'%vxi)
	ctrl.write('set %s.1.20 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
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

def stock_add_to_all(cd='1c', ce='2c', he='2h'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(cd, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(ce, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(he, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd, lac(p['CA'][0], p['CA'][0], i), lac(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ce, lac(p['CA'][0], p['0_C'][0], i), lac(p['CA'][1], p['0_C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he, lac(p['HA'][0], p['0_H'][0], i), lac(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd, 'HA'), lac(p['CA_HA'][0], p['CA_HA'][0], i), lac(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CA', cd), lac(p['CA_CA'][0], p['CA_CA'][0], i), lac(p['CA_CA'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce, he), lac(p['CA_HA'][0], p['HA_mCA'][0], i), lac(p['CA_HA'][1], p['HA_mCA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce, cd), lac(p['CA_CA'][0], p['CA_mCA'][0], i), lac(p['CA_CA'][1], p['CA_mCA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CA', ce), lac(p['CA_CA'][0], p['CA_mCA'][0], i), lac(p['CA_CA'][1], p['CA_mCA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he, ce, cd), lac(p['F_CA_CA_HA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_HA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HA', cd, ce), lac(p['F_CA_CA_HA'][0], p['A_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['A_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HA', 'CA', ce), lac(p['F_CA_CA_HA'][0], p['A_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['A_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', ce, he), lac(p['F_CA_CA_HA'][0], p['CA_Close'][0], i), lac(p['F_CA_CA_HA'][1], p['CA_Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', ce, cd), lac(p['F_CA_CA_CA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_CA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce, 'CA', ce), lac(p['F_CA_CA_CA'][0], p['A_CA_CA_CA'][0], i), lac(p['F_CA_CA_CA'][1], p['A_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd, 'CA', cd), lac(p['F_CA_CA_CA'][0], p['A_CA_CA_CA'][0], i), lac(p['F_CA_CA_CA'][1], p['A_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HA', cd, ce), lac(p['F_CA_CA_HA'][0], p['A_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['A_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', cd, ce), lac(p['F_CA_CA_CA'][0], p['A_CA_CA_CA'][0], i), lac(p['F_CA_CA_CA'][1], p['A_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', cd, 'HA'), lac(p['F_CA_CA_HA'][0], p['A_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['A_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CA', cd), lac(p['F_CT_CA_CA'][0], p['A_CT_CA_CA'][0], i), lac(p['F_CT_CA_CA'][1], p['A_CT_CA_CA'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ce, cd, 'CA'), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ce, cd, 'HA'), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ce, 'CA', 'HA'), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, ce, 'CA', 'HA'), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', cd, ce, 'HA'), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ce, 'CA', ce), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce, cd, 'CA', cd), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', ce, cd, 'HA'), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce, 'CA', ce, cd), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HA', cd, 'CA', cd), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', cd, ce, 'CA'), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CA', cd, ce), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CA', cd, 'HA'), lac(p['Ring_Dihe'][0], p['Ring_Dihe'][0], i), lac(p['Ring_Dihe'][1], p['Ring_Dihe'][1], i), lac(p['Ring_Dihe'][2], p['Ring_Dihe'][2], i), lac(p['Ring_Dihe'][3], p['Ring_Dihe'][3], i))
	#	Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', ce, he), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
	#	Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, 'CA', cd, 'CT'), lac(p['Ring_imp'][0], p['Ring_imp'][0], i), lac(p['Ring_imp'][1], p['Ring_imp'][1], i), lac(p['Ring_imp'][2], p['Ring_imp'][2], i))
#		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', cd, 'HA'), lac(p['Ring_imp'][0], p['Ring_imp'][0], i), lac(p['Ring_imp'][1], p['Ring_imp'][1], i), lac(p['Ring_imp'][2], p['Ring_imp'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd, lac(p['CA'][2], p['CA'][2], i), lac(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ce, lac(p['CA'][2], p['0_C'][2], i), lac(p['CA'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he, lac(p['HA'][2], p['0_H'][2], i), lac(p['HA'][3], p['0_H'][3], i))
