# PHE to NAL Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/PHE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/NAL.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                changeLJPair(parm, ':{}@HH2 :{}@HZ2 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HH2 :{}@CE2 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HT :{}@HH1 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HT :{}@CZ1 0 0'.format(vxi, vxi)).execute()
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
                change(parm, 'charge', ':{}@CE1'.format(vxi), bc['CE1']+((fc['CE1']-bc['CE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), bc['HE1']+((fc['HE1']-bc['HE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ1'.format(vxi), bc['CZ']+((fc['CE2']-bc['CZ'])/10)*i).execute()
                change(parm, 'charge', ':{}@CH1'.format(vxi), bc['HZ']+((fc['CZ2']-bc['HZ'])/10)*i).execute()
                change(parm, 'charge', ':{}@HH1'.format(vxi), (fc['HH1']/10)*i).execute()
                change(parm, 'charge', ':{}@CT'.format(vxi), (fc['CT']/10)*i).execute()
                change(parm, 'charge', ':{}@HT'.format(vxi), (fc['HT']/10)*i).execute()
                change(parm, 'charge', ':{}@CH2'.format(vxi), (fc['CH2']/10)*i).execute()
                change(parm, 'charge', ':{}@HH2'.format(vxi), (fc['HH2']/10)*i).execute()
                change(parm, 'charge', ':{}@CZ2'.format(vxi), bc['HE2']+((fc['CZ2']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vxi), (fc['HZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@CE2'.format(vxi), bc['CE2']+((fc['CE2']-bc['CE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']+((fc['CD2']-bc['CD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	CZ = struct.residue_dict[aa].atom_dict['CZ']
	HZ = struct.residue_dict[aa].atom_dict['HZ']
	CE2 = struct.residue_dict[aa].atom_dict['CE2']
	HE2 = struct.residue_dict[aa].atom_dict['HE2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'HE2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ2'))
				pdb.write(atom.superimposed1('HZ2', CE2))
			elif atom.get_name() == 'CZ' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ1'))
			elif atom.get_name() == 'HZ' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CH1'))
				pdb.write(atom.superimposed1('HH1', CZ))
				pdb.write(atom.thirdone_between('CT', HZ, HE2))
				pdb.write(atom.superimposed1('HT', HZ))
				pdb.write(atom.thirdtwo_between('CH2', HZ, HE2))
				pdb.write(atom.superimposed1('HH2', HE2))
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
	ch1 = var[0]
	hh1 = var[1]
	ct = var[2]
	ht = var[3]
	ch2 = var[4]
	hh2 = var[5]
	cz2 = var[6]
	hz2 = var[7]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/PHE-NAL.pdb\n"%vxi)
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
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "C"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "C"\n'%vxi)
	ctrl.write('set %s.1.23 element "C"\n'%vxi)
	ctrl.write('set %s.1.24 element "H"\n'%vxi)
	ctrl.write('set %s.1.25 element "C"\n'%vxi)
	ctrl.write('set %s.1.26 element "O"\n'%vxi)
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
	ctrl.write('set %s.1.13 name "CZ1"\n'%vxi)
	ctrl.write('set %s.1.14 name "CH1"\n'%vxi)
	ctrl.write('set %s.1.15 name "HH1"\n'%vxi)
	ctrl.write('set %s.1.16 name "CT"\n'%vxi)
	ctrl.write('set %s.1.17 name "HT"\n'%vxi)
	ctrl.write('set %s.1.18 name "CH2"\n'%vxi)
	ctrl.write('set %s.1.19 name "HH2"\n'%vxi)
	ctrl.write('set %s.1.20 name "CZ2"\n'%vxi)
	ctrl.write('set %s.1.21 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.22 name "CE2"\n'%vxi)
	ctrl.write('set %s.1.23 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.24 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.25 name "C"\n'%vxi)
	ctrl.write('set %s.1.26 name "O"\n'%vxi)
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
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, ch1))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, hh1))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, ct))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, ht))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, ch2))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, hh2))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, cz2))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, hz2))
	ctrl.write('set %s.1.22 type "CA"\n'%vxi)
	ctrl.write('set %s.1.23 type "CA"\n'%vxi)
	ctrl.write('set %s.1.24 type "HA"\n'%vxi)
	ctrl.write('set %s.1.25 type "C"\n'%vxi)
	ctrl.write('set %s.1.26 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.23\n'%(vxi, vxi))
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
	ctrl.write('bond %s.1.22 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.26\n'%(vxi, vxi))
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
	ch1 = var[0]
	hh1 = var[1]
	ct = var[2]
	ht = var[3]
	ch2 = var[4]
	hh2 = var[5]
	cz2 = var[6]
	hz2 = var[7]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(ch1, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hh1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ct, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(ht, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ch2, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hh2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cz2, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hz2, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ch1, cal(p['HA'][0], p['CA'][0], i), cal(p['HA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh1, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ct, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ht, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ch2, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh2, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz2, cal(p['HA'][0], p['CA'][0], i), cal(p['HA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz2, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CA', ch1), cal(p['CA_HA'][0], p['CA_CA'][0], i), cal(p['CA_HA'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch1, hh1), cal(p['CA_HA'][0], p['CA_HA'][0], i), cal(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch1, ct), cal(p['CA_tCA4'][0], p['CA_CA'][0], i), cal(p['CA_tCA4'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ct, ht), cal(p['HA_tCA4'][0], p['CA_HA'][0], i), cal(p['HA_tCA4'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ct, ch2), cal(p['CA_tCA4'][0], p['CA_CA'][0], i), cal(p['CA_tCA4'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch2, hh2), cal(p['HA_tCA4'][0], p['CA_HA'][0], i), cal(p['HA_tCA4'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch2, cz2), cal(p['CA_tCA4'][0], p['CA_CA'][0], i), cal(p['CA_tCA4'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz2, 'CA'), cal(p['CA_HA'][0], p['CA_CA'][0], i), cal(p['CA_HA'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz2, hz2), cal(p['CA_HA'][0], p['CA_HA'][0], i), cal(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', 'CA', ch1), cal(p['F_CA_CA_HA'][0], p['F_CA_CA_CA'][0], i), cal(p['F_CA_CA_HA'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', ch1, hh1), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', ch1, ct), cal(p['A_Phe'][0], p['F_CA_CA_CA'][0], i), cal(p['A_Phe'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh1, ch1, ct), cal(p['A_Phe'][0], p['F_CA_CA_HA'][0], i), cal(p['A_Phe'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch1, ct, ht), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch1, ct, ch2), cal(p['Dritt'][0], p['F_CA_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ht, ct, ch2), cal(p['Dritt'][0], p['F_CA_CA_HA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ct, ch2, hh2), cal(p['Dritt'][0], p['F_CA_CA_HA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ct, ch2, cz2), cal(p['Dritt'][0], p['F_CA_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh2, ch2, cz2), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch2, cz2, 'CA'), cal(p['A_Phe'][0], p['F_CA_CA_CA'][0], i), cal(p['A_Phe'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch2, cz2, hz2), cal(p['A_Phe'][0], p['F_CA_CA_HA'][0], i), cal(p['A_Phe'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz2, cz2, 'CA'), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz2, 'CA', 'CA'), cal(p['F_CA_CA_HA'][0], p['F_CA_CA_CA'][0], i), cal(p['F_CA_CA_HA'][1], p['F_CA_CA_CA'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', ch1, hh1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', ch1, ct), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', ch1, ct, ht), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', ch1, ct, ch2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh1, ch1, ct, ht), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh1, ch1, ct, ch2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch1, ct, ch2, hh2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch1, ct, ch2, cz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ht, ct, ch2, hh2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ht, ct, ch2, cz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ct, ch2, cz2, hz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ct, ch2, cz2, 'CA'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh2, ch2, cz2, hz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh2, ch2, cz2, 'CA'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch2, cz2, 'CA', 'CA'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz2, cz2, 'CA', 'CA'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_0'][2], i), cal(p['Ring_0'][3], p['Ring_0'][3], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', ch1, hh1), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', ct, ht), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', ch2, hh2), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', cz2, hz2), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ch1, cal(p['HA'][2], p['CA'][2], i), cal(p['HA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh1, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ct, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ht, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ch2, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh2, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz2, cal(p['HA'][2], p['CA'][2], i), cal(p['HA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz2, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
