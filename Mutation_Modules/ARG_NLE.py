# NLE to ARG Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/ARG.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/NLE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HE1'.format(vxi), ':{}@HH11'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@HE1'.format(vxi), ':{}@HH21'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@HE2'.format(vxi), ':{}@HH11'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@HE2'.format(vxi), ':{}@HH21'.format(vxi), '0', '0').execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), bc['CD']+((fc['CD']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), bc['HD3']+((fc['HD3']-bc['HD3'])/10)*i).execute()
                change(parm, 'charge', ':{}@NE'.format(vxi), bc['NE']-(bc['NE']/10)*i).execute()
                change(parm, 'charge', ':{}@HE'.format(vxi), bc['HE']-(bc['HE']/10)*i).execute()
                change(parm, 'charge', ':{}@CE'.format(vxi), bc['CZ']+((fc['CE']-bc['CZ'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), (fc['HE1']/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), (fc['HE2']/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), (fc['HE3']/10)*i).execute()
                change(parm, 'charge', ':{}@NH1'.format(vxi), bc['NH1']-(bc['NH1']/10)*i).execute()
                change(parm, 'charge', ':{}@HH11'.format(vxi), bc['HH11']-(bc['HH11']/10)*i).execute()
                change(parm, 'charge', ':{}@HH12'.format(vxi), bc['HH12']-(bc['HH12']/10)*i).execute()
                change(parm, 'charge', ':{}@NH2'.format(vxi), bc['NH2']-(bc['NH2']/10)*i).execute()
                change(parm, 'charge', ':{}@HH21'.format(vxi), bc['HH21']-(bc['HH21']/10)*i).execute()
                change(parm, 'charge', ':{}@HH22'.format(vxi), bc['HH22']-(bc['HH22']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
                d = netCharge(parm).execute()
                change(parm, 'charge', ':PN', '{:.3f}'.format(-d)).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        NE = struct.residue_dict[aa].atom_dict['NE']
        NH1 = struct.residue_dict[aa].atom_dict['NH1']
        NH2 = struct.residue_dict[aa].atom_dict['NH2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'CZ' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CE')) 
				pdb.write(atom.superimposed1('HE1', NH1)) 
				pdb.write(atom.superimposed1('HE2', NH2)) 
				pdb.write(atom.superimposed1('HE3', NE)) 
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
	hd1 = var[0]
	ne = var[1]
	he = var[2]
	cz = var[3]
	he1 = var[4]
	he2 = var[5]
	he3 = var[6]
	nh1 = var[7]
	hh1 = var[8]
	nh2 = var[9]
	hh2 = var[10]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/NLE-ARG.pdb\n"%vxi)
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
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "N"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "N"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "N"\n'%vxi)
	ctrl.write('set %s.1.24 element "H"\n'%vxi)
	ctrl.write('set %s.1.25 element "H"\n'%vxi)
	ctrl.write('set %s.1.26 element "C"\n'%vxi)
	ctrl.write('set %s.1.27 element "O"\n'%vxi)
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
	ctrl.write('set %s.1.11 name "CD"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.13 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.14 name "NE"\n'%vxi)
	ctrl.write('set %s.1.15 name "HE"\n'%vxi)
	ctrl.write('set %s.1.16 name "CE"\n'%vxi)
	ctrl.write('set %s.1.17 name "HE1"\n'%vxi)
	ctrl.write('set %s.1.18 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.19 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.20 name "NH1"\n'%vxi)
	ctrl.write('set %s.1.21 name "HH11"\n'%vxi)
	ctrl.write('set %s.1.22 name "HH12"\n'%vxi)
	ctrl.write('set %s.1.23 name "NH2"\n'%vxi)
	ctrl.write('set %s.1.24 name "HH21"\n'%vxi)
	ctrl.write('set %s.1.25 name "HH22"\n'%vxi)
	ctrl.write('set %s.1.26 name "C"\n'%vxi)
	ctrl.write('set %s.1.27 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "CT"\n'%vxi)
	ctrl.write('set %s.1.9 type "HC"\n'%vxi)
	ctrl.write('set %s.1.10 type "HC"\n'%vxi)
	ctrl.write('set %s.1.11 type "CT"\n'%vxi)
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, hd1))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, hd1))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, ne))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, he))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, cz))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, he1))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, he2))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, he3))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, nh1))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, hh1))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, hh1))
	ctrl.write('set %s.1.23 type "%s"\n'%(vxi, nh2))
	ctrl.write('set %s.1.24 type "%s"\n'%(vxi, hh2))
	ctrl.write('set %s.1.25 type "%s"\n'%(vxi, hh2))
	ctrl.write('set %s.1.26 type "C"\n'%vxi)
	ctrl.write('set %s.1.27 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.26\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.26 %s.1.27\n'%(vxi, vxi))
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
	hd1 = var[0]
	ne = var[1]
	he = var[2]
	cz = var[3]
	he1 = var[4]
	he2 = var[5]
	he3 = var[6]
	nh1 = var[7]
	hh1 = var[8]
	nh2 = var[9]
	hh2 = var[10]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(hd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ne, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(he, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cz, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(he1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(he2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(he3, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(nh1, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(hh1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(nh2, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(hh2, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd1, lac(p['HC'][0], p['H1'][0], i), lac(p['HC'][1], p['H1'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ne, lac(p['0_N'][0], p['NA'][0], i), lac(p['0_N'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz, lac(p['CT'][0], p['CA'][0], i), lac(p['CT'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he1, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he2, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he3, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), nh1, lac(p['0_N'][0], p['NA'][0], i), lac(p['0_N'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh1, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), nh2, lac(p['0_N'][0], p['NA'][0], i), lac(p['0_N'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh2, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hd1), lac(p['CT_HC'][0], p['CT_HC'][0], i), lac(p['CT_HC'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', ne), lac(p['N2_mH'][0], p['CT_N'][0], i), lac(p['N2_mH'][1], p['CT_N'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne, he), lac(p['H_mHC'][0], p['NA_H'][0], i), lac(p['H_mHC'][1], p['NA_H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne, cz), lac(p['N2_mH2'][0], p['CA_N'][0], i), lac(p['N2_mH2'][1], p['CA_N'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, he1), lac(p['CT_HC'][0], p['HC_sN2'][0], i), lac(p['CT_HC'][1], p['HC_sN2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, he2), lac(p['CT_HC'][0], p['HC_sN2'][0], i), lac(p['CT_HC'][1], p['HC_sN2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, he3), lac(p['CT_HC'][0], p['HC_sN2'][0], i), lac(p['CT_HC'][1], p['HC_sN2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, nh1), lac(p['N2_mH2'][0], p['CA_N'][0], i), lac(p['N2_mH2'][1], p['CA_N'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, nh2), lac(p['N2_mH2'][0], p['CA_N'][0], i), lac(p['N2_mH2'][1], p['CA_N'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(nh1, hh1), lac(p['H_mHC'][0], p['NA_H'][0], i), lac(p['H_mHC'][1], p['NA_H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(nh2, hh2), lac(p['H_mHC'][0], p['NA_H'][0], i), lac(p['H_mHC'][1], p['NA_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hd1), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, 'CT', hd1), lac(p['H_C_H'][0], p['H_C_H'][0], i), lac(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', ne), lac(p['C_C_C'][0], p['C_C_N2'][0], i), lac(p['C_C_C'][1], p['C_C_N2'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, 'CT', ne), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ne, he), lac(p['Dritt'][0], p['C_N_H'][0], i), lac(p['Dritt'][1], p['C_N_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ne, cz), lac(p['Dritt'][0], p['C_N_C'][0], i), lac(p['Dritt'][1], p['C_N_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he, ne, cz), lac(p['Close'][0], p['F_CA_CA_HA'][0], i), lac(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, he1), lac(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), lac(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, he2), lac(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), lac(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, he3), lac(p['C_C_H'][0], p['Close'][0], i), lac(p['C_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, nh1), lac(p['C_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['C_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, nh2), lac(p['C_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['C_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he1, cz, nh1), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he1, cz, nh2), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he2, cz, nh1), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he2, cz, nh2), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he3, cz, nh1), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he3, cz, nh2), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(nh1, cz, nh2), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz, nh1, hh1), lac(p['Dritt'][0], p['F_CA_CA_HA'][0], i), lac(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz, nh2, hh2), lac(p['Dritt'][0], p['F_CA_CA_HA'][0], i), lac(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh1, nh1, hh1), lac(p['Close'][0], p['H_N_H'][0], i), lac(p['Close'][1], p['H_N_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh2, nh2, hh2), lac(p['Close'][0], p['H_N_H'][0], i), lac(p['Close'][1], p['H_N_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he1, cz, he2), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he1, cz, he3), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he2, cz, he3), lac(p['H_C_H'][0], p['F_CA_CA_CA'][0], i), lac(p['H_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ne, he), lac(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), lac(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), lac(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), lac(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, 'CT', ne, he), lac(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), lac(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), lac(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), lac(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ne, cz), lac(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), lac(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), lac(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), lac(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, 'CT', ne, cz), lac(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), lac(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), lac(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), lac(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, he1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, he2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, he3), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ne, cz, he1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ne, cz, he2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ne, cz, he3), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, nh1), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, nh2), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ne, cz, nh1), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he, ne, cz, nh2), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne, cz, nh1, hh1), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne, cz, nh2, hh2), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he1, cz, nh1, hh1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he1, cz, nh2, hh2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he2, cz, nh1, hh1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he2, cz, nh2, hh2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, cz, nh1, hh1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, cz, nh2, hh2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nh1, cz, nh2, hh2), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nh2, cz, nh1, hh1), lac(p['Ring_0'][0], p['X_CA_N_X'][0], i), lac(p['Ring_0'][1], p['X_CA_N_X'][1], i), lac(p['Ring_0'][2], p['X_CA_N_X'][2], i), lac(p['Ring_0'][3], p['X_CA_N_X'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz, hh1, nh1, hh1), lac(p['Imp_0'][0], p['Ami_imp'][0], i), lac(p['Imp_0'][1], p['Ami_imp'][1], i), lac(p['Imp_0'][2], p['Ami_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz, hh2, nh2, hh2), lac(p['Imp_0'][0], p['Ami_imp'][0], i), lac(p['Imp_0'][1], p['Ami_imp'][1], i), lac(p['Imp_0'][2], p['Ami_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne, nh1, cz, nh2), lac(p['Imp_0'][0], p['Car_imp'][0], i), lac(p['Imp_0'][1], p['Car_imp'][1], i), lac(p['Imp_0'][2], p['Car_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, he), lac(p['Imp_0'][0], p['Car_imp'][0], i), lac(p['Imp_0'][1], p['Car_imp'][1], i), lac(p['Imp_0'][2], p['Car_imp'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd1, lac(p['HC'][2], p['H1'][2], i), lac(p['HC'][3], p['H1'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ne, lac(p['0_N'][2], p['NA'][2], i), lac(p['0_N'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz, lac(p['CT'][2], p['CA'][2], i), lac(p['CT'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he1, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he2, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he3, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), nh1, lac(p['0_N'][2], p['NA'][2], i), lac(p['0_N'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh1, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), nh2, lac(p['0_N'][2], p['NA'][2], i), lac(p['0_N'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh2, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
