# ANN to ILE Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/ANN.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/ILE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HB3'.format(vxi), ':{}@HG21'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@HD1'.format(vxi), ':{}@HD2'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@HG12'.format(vxi), ':{}@HD12'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@HG13'.format(vxi), ':{}@HD13'.format(vxi), '0', '0').execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB'.format(vxi), bc['HB2']+((fc['HB']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']-(bc['HB3']/10)*i).execute()
                change(parm, 'charge', ':{}@CG2'.format(vxi), fc['CG2']/10*i).execute()
                change(parm, 'charge', ':{}@HG21'.format(vxi), fc['HG21']/10*i).execute()
                change(parm, 'charge', ':{}@HG22'.format(vxi), fc['HG22']/10*i).execute()
                change(parm, 'charge', ':{}@HG23'.format(vxi), fc['HG23']/10*i).execute()
                change(parm, 'charge', ':{}@CG1'.format(vxi), bc['CG']+((fc['CG1']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG12'.format(vxi), (fc['HG12']/10)*i).execute()
                change(parm, 'charge', ':{}@HG13'.format(vxi), (fc['HG13']/10)*i).execute()
                change(parm, 'charge', ':{}@CD3'.format(vxi), bc['CD1']-(bc['CD1']/10)*i).execute()
                change(parm, 'charge', ':{}@HD1'.format(vxi), bc['HD1']-(bc['HD1']/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']-(bc['CD2']/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), bc['HD2']-(bc['HD2']/10)*i).execute()
                change(parm, 'charge', ':{}@CD1'.format(vxi), bc['CE']+((fc['CD1']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD11'.format(vxi), bc['HE']+((fc['HD11']-bc['HE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD12'.format(vxi), (fc['HD12']/10)*i).execute()
                change(parm, 'charge', ':{}@HD13'.format(vxi), (fc['HD13']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CB = struct.residue_dict[aa].atom_dict['CB']
        HB3 = struct.residue_dict[aa].atom_dict['HB3']
        CD1 = struct.residue_dict[aa].atom_dict['CD1']
        CD2 = struct.residue_dict[aa].atom_dict['CD2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'HB2' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HB'))
			elif atom.get_name() == 'HB3' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.halfway_between1('CG2', CB, HB3)) 
				pdb.write(atom.superimposed1('HG21', HB3)) 
				pdb.write(atom.superimposed2('HG22', HB3)) 
				pdb.write(atom.superimposed3('HG23', HB3)) 
			elif atom.get_name() == 'CG' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CG1'))
				pdb.write(atom.superimposed1('HG12', CD1)) 
				pdb.write(atom.superimposed1('HG13', CD2)) 
			elif atom.get_name() == 'CD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CD3')) 
			elif atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CD1')) 
			elif atom.get_name() == 'HE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HD11')) 
				pdb.write(atom.superimposed2('HD12', CD1)) 
				pdb.write(atom.superimposed2('HD13', CD2)) 
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
	return var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14

def lib_make(ff, outputfile, vxi='VXI', var=variablemake()):
	metcar = var[0]
	methyd = var[1]
	hydhyd = var[2]
	cg = var[3]
	hg1 = var[4]
	hg2 = var[5]
	cd1 = var[6]
	hd11 = var[7]
	cd2 = var[8]
	hd21 = var[9]
	cd = var[10]
	hd1 = var[11]
	hd2 = var[12]
	hd3 = var[13]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ILE-ANN.pdb\n"%vxi)
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
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "C"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "C"\n'%vxi)
	ctrl.write('set %s.1.24 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG2"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG21"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG22"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG23"\n'%vxi)
	ctrl.write('set %s.1.12 name "CG1"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG12"\n'%vxi)
	ctrl.write('set %s.1.14 name "HG13"\n'%vxi)
	ctrl.write('set %s.1.15 name "CD3"\n'%vxi)
	ctrl.write('set %s.1.16 name "HD1"\n'%vxi)
	ctrl.write('set %s.1.17 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.18 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.19 name "CD1"\n'%vxi)
	ctrl.write('set %s.1.20 name "HD11"\n'%vxi)
	ctrl.write('set %s.1.21 name "HD12"\n'%vxi)
	ctrl.write('set %s.1.22 name "HD13"\n'%vxi)
	ctrl.write('set %s.1.23 name "C"\n'%vxi)
	ctrl.write('set %s.1.24 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, hydhyd))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, cg))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, hg1))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, hg2))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, cd1))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, hd11))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, cd2))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, hd21))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, cd))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, hd1))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, hd2))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, hd3))
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
	ctrl.write('bond %s.1.5 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.22\n'%(vxi, vxi))
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
	metcar = var[0]
	methyd = var[1]
	hydhyd = var[2]
	cg = var[3]
	hg1 = var[4]
	hg2 = var[5]
	cd1 = var[6]
	hd11 = var[7]
	cd2 = var[8]
	hd21 = var[9]
	cd = var[10]
	hd1 = var[11]
	hd2 = var[12]
	hd3 = var[13]
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
        Frcmod_creator.TYPE_insert(metcar, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(methyd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(hydhyd, 'H', 'sp3')
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
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', cg), cal(p['CT_CA'][0], p['CT_CT'][0], i), cal(p['CT_CA'][1], p['CT_CT'][1], i))
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
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(metcar, 'CT', cg), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', cg), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
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
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', cg, hg1), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', cg, hg2), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar, 'CT', cg, hg1), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar, 'CT', cg, hg2), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar, 'CT', cg, cd1), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar, 'CT', cg, cd2), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(methyd, metcar, 'CT', cg), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, hg1), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, hg2), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd1), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd2), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', cg, cd1), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', cg, cd2), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd1), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd2), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, cd1), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, cd2), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd1, hd11), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, hd21), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd1, cd), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, cd), cal(p['Ring_Dihe'][0], p['Ring_0'][0], i), cal(p['Ring_Dihe'][1], p['Ring_0'][1], i), cal(p['Ring_Dihe'][2], p['Ring_0'][2], i), cal(p['Ring_Dihe'][3], p['Ring_0'][3], i))
		#Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X', 'X', cd, hd1), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
		#Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X', 'X', cd1, hd11), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
		#Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X', 'X', cd2, hd21), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
		#Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, cg, cd1, 'CT'), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
		#Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X', 'X', cd, 'HA'), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
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

                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][0], p['0_C'][0], i), lac(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][0], p['HC'][0], i), lac(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), lac(p['CT_CT'][0], p['CT_mH'][0], i), lac(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), lac(p['HC_sC'][0], p['CT_HC'][0], i), lac(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), lac(p['CT_HC'][0], p['HC_mH'][0], i), lac(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), lac(p['C_C_H'][0], p['Dritt'][0], i), lac(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), lac(p['H_C_H'][0], p['Close'][0], i), lac(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), lac(p['C_C_C'][0], p['C_C_C'][0], i), lac(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', metcar), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', hydhyd), lac(p['H_C_H'][0], p['H_C_H'][0], i), lac(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', metcar), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), lac(p['C_C_C_H'][0], p['0_1'][0], i), lac(p['C_C_C_H'][1], p['0_1'][1], i), lac(p['C_C_C_H'][2], p['0_1'][2], i), lac(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', metcar, methyd), lac(p['H_C_C_H'][0], p['0_1'][0], i), lac(p['H_C_C_H'][1], p['0_1'][1], i), lac(p['H_C_C_H'][2], p['0_1'][2], i), lac(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', metcar, methyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][2], p['0_C'][2], i), lac(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][2], p['HC'][2], i), lac(p['0_H'][3], p['HC'][3], i))
