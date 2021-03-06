# ILE to VAL Mutation

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
        with open('Param_files/AminoAcid/NVA.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HD11 :{}@HD21 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HG2 :{}@HD2 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HG3 :{}@HD3 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), (fc['HG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), (fc['HG3']/10)*i).execute()
                change(parm, 'charge', ':{}@CD1'.format(vxi), bc['CD1']-(bc['CD1']/10)*i).execute()
                change(parm, 'charge', ':{}@HD11'.format(vxi), bc['HD1']-(bc['HD1']/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']-(bc['CD2']/10)*i).execute()
                change(parm, 'charge', ':{}@HD21'.format(vxi), bc['HD2']-(bc['HD2']/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), bc['CE']+((fc['CD']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD1'.format(vxi), bc['HE']+((fc['HD1']-bc['HE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), (fc['HD2']/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), (fc['HD3']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CD1 = struct.residue_dict[aa].atom_dict['CD1']
        CD2 = struct.residue_dict[aa].atom_dict['CD2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'CG' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.superimposed1('HG2', CD1)) 
				pdb.write(atom.superimposed1('HG3', CD2)) 
			elif atom.get_name() == 'HD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HD11')) 
			elif atom.get_name() == 'HD2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HD21')) 
			elif atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CD')) 
			elif atom.get_name() == 'HE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HD1')) 
				pdb.write(atom.superimposed2('HD2', CD1)) 
				pdb.write(atom.superimposed2('HD3', CD2)) 
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

def lib_make(ff, outputfile, vxi='VXI', cg='1c', hg1='0h', hg2='1h', cd1='2c', hd11='2h', cd2='3c', hd21='3h', cd='4c', hd1='4h', hd2='5h', hd3='6h'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ANN-NVA.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
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
	ctrl.write('set %s.1.9 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.11 name "CD1"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD11"\n'%vxi)
	ctrl.write('set %s.1.13 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.14 name "HD21"\n'%vxi)
	ctrl.write('set %s.1.15 name "CD"\n'%vxi)
	ctrl.write('set %s.1.16 name "HD1"\n'%vxi)
	ctrl.write('set %s.1.17 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.18 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.19 name "C"\n'%vxi)
	ctrl.write('set %s.1.20 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, cg))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, hg1))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, hg2))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, cd1))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, hd11))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, cd2))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, hd21))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, cd))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, hd1))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, hd2))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, hd3))
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
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.18\n'%(vxi, vxi))
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

def stock_add_to_all(cg='1c', hg1='0h', hg2='1h', cd1='2c', hd11='2h', cd2='3c', hd21='3h', cd='4c', hd1='4h', hd2='5h', hd3='6h'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(cg, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hg1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hg2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cd1, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hd11, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cd2, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hd21, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cd, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hd2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hd3, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cg, cal(p['CA'][0], p['CT'][0], i), cal(p['CA'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hg1, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hg2, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd1, cal(p['CA'][0], p['0_C'][0], i), cal(p['CA'][1], p['0_C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd11, cal(p['HA'][0], p['0_H'][0], i), cal(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd2, cal(p['CA'][0], p['0_C'][0], i), cal(p['CA'][1], p['0_C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd21, cal(p['HA'][0], p['0_H'][0], i), cal(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd, cal(p['CA'][0], p['CT'][0], i), cal(p['CA'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd1, cal(p['HA'][0], p['HC'][0], i), cal(p['HA'][1], p['HC'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd2, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd3, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', cg), cal(p['CT_CA'][0], p['CT_CT'][0], i), cal(p['CA_CA'][1], p['CT_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cg, hg1), cal(p['HC_sCA'][0], p['CT_HC'][0], i), cal(p['HC_sCA'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cg, hg2), cal(p['HC_sCA'][0], p['CT_HC'][0], i), cal(p['HC_sCA'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd, hd1), cal(p['CA_HA'][0], p['CT_HC'][0], i), cal(p['CA_HA'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd, hd2), cal(p['HC_sCA'][0], p['CT_HC'][0], i), cal(p['HC_sCA'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd, hd3), cal(p['HC_sCA'][0], p['CT_HC'][0], i), cal(p['HC_sCA'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd1, hd11), cal(p['CA_HA'][0], p['HA_mCT'][0], i), cal(p['CA_HA'][1], p['HA_mCT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd2, hd21), cal(p['CA_HA'][0], p['HA_mCT'][0], i), cal(p['CA_HA'][1], p['HA_mCT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cg, cd1), cal(p['CA_CA'][0], p['CA_mCT'][0], i), cal(p['CA_CA'][1], p['CA_mCT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cg, cd2), cal(p['CA_CA'][0], p['CA_mCT'][0], i), cal(p['CA_CA'][1], p['CA_mCT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd1, cd), cal(p['CA_CA'][0], p['CA_mCT'][0], i), cal(p['CA_CA'][1], p['CA_mCT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd2, cd), cal(p['CA_CA'][0], p['CA_mCT'][0], i), cal(p['CA_CA'][1], p['CA_mCT'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, cd, hd2), cal(p['A_CA_CA_HA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, cd, hd3), cal(p['A_CA_CA_HA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd2, cd, hd3), cal(p['A_CA_CA_CA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_CA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, cd, cd1), cal(p['A_CA_CA_HA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, cd, cd2), cal(p['A_CA_CA_HA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd11, cd1, cd), cal(p['A_CA_CA_HA'][0], p['Close'][0], i), cal(p['A_CA_CA_HA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd21, cd2, cd), cal(p['A_CA_CA_HA'][0], p['Close'][0], i), cal(p['A_CA_CA_HA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd11, cd1, cg), cal(p['A_CA_CA_HA'][0], p['Dritt'][0], i), cal(p['A_CA_CA_HA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd21, cd2, cg), cal(p['A_CA_CA_HA'][0], p['Dritt'][0], i), cal(p['A_CA_CA_HA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd2, cd, cd1), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd2, cd, cd2), cal(p['A_CA_CA_CA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_CA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd3, cd, cd1), cal(p['A_CA_CA_CA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_CA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd3, cd, cd2), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd1, cg, cd2), cal(p['A_CA_CA_CA'][0], p['Close'][0], i), cal(p['A_CA_CA_CA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd, cd1, cg), cal(p['A_CA_CA_CA'][0], p['Dritt'][0], i), cal(p['A_CA_CA_CA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd1, cd, cd2), cal(p['A_CA_CA_CA'][0], p['Close'][0], i), cal(p['A_CA_CA_CA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd, cd2, cg), cal(p['A_CA_CA_CA'][0], p['Dritt'][0], i), cal(p['A_CA_CA_CA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hg1, cg, cd1), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hg1, cg, cd2), cal(p['A_CA_CA_CA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_CA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hg2, cg, cd1), cal(p['A_CA_CA_CA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_CA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hg2, cg, cd2), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hg1, cg, hg2), cal(p['A_CA_CA_CA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_CA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', cg), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hg1, cg, 'CT'), cal(p['A_CA_CA_HA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hg2, cg, 'CT'), cal(p['A_CA_CA_HA'][0], p['C_C_H'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd1, cg, 'CT'), cal(p['A_CA_CA_HA'][0], p['C_C_C'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd2, cg, 'CT'), cal(p['A_CA_CA_HA'][0], p['C_C_C'][0], i), cal(p['A_CA_CA_HA'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', cg), cal(p['C_C_CA'][0], p['C_C_H'][0], i), cal(p['C_C_CA'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, cd, cd1, hd11), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd2, cd, cd1, hd11), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd3, cd, cd1, hd11), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, cd, cd2, hd21), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd2, cd, cd2, hd21), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd3, cd, cd2, hd21), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd1, cd, cd2, hd21), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, cd, cd1, hd11), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd1, cg, cd2, hd21), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, cg, cd1, hd11), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg1, cg, cd1, hd11), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg1, cg, cd2, hd21), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg2, cg, cd1, hd11), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg2, cg, cd2, hd21), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg1, cg, cd1, cd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg1, cg, cd2, cd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg2, cg, cd1, cd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hg2, cg, cd2, cd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd1, cd, hd1), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd1, cd, hd2), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd1, cd, hd3), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, cd, hd1), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, cd, hd2), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, cd, hd3), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, cd, cd1), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd1, cd, cd2), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd1, cg, cd2, cd), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, cg, cd1, cd), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, hg1), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, hg2), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, hg1), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, hg2), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd1), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd2), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd1), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd2), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, cd1), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, cd2), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd1, hd11), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, hd21), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd1, cd), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, cd), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
#         #       Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X', 'X', ce, he), cal(p['Ring_imp'][0], p['Ring_imp_0'][0], i), cal(p['Ring_imp'][1], p['Ring_imp_0'][1], i), cal(p['Ring_imp'][2], p['Ring_imp_0'][2], i))
#               Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, 'CA', cd, 'CT'), cal(p['Ring_imp'][0], p['Ring_imp'][0], i), cal(p['Ring_imp'][1], p['Ring_imp'][1], i), cal(p['Ring_imp'][2], p['Ring_imp'][2], i))
#               Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X', 'X', cd, 'HA'), cal(p['Ring_imp'][0], p['Ring_imp'][0], i), cal(p['Ring_imp'][1], p['Ring_imp'][1], i), cal(p['Ring_imp'][2], p['Ring_imp'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cg, cal(p['CA'][2], p['CT'][2], i), cal(p['CA'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hg1, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hg2, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd1, cal(p['CA'][2], p['0_C'][2], i), cal(p['CA'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd11, cal(p['HA'][2], p['0_H'][2], i), cal(p['HA'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd2, cal(p['CA'][2], p['0_C'][2], i), cal(p['CA'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd21, cal(p['HA'][2], p['0_H'][2], i), cal(p['HA'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd, cal(p['CA'][2], p['CT'][2], i), cal(p['CA'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd1, cal(p['HA'][2], p['HC'][2], i), cal(p['HA'][3], p['HC'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd2, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd3, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
