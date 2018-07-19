# LYS to ARG Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/LYS.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/ARG.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HZ1 :{}@HH11 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HZ1 :{}@HH21 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HZ2 :{}@HH11 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HZ2 :{}@HH21 0 0'.format(vxi, vxi)).execute()
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
                change(parm, 'charge', ':{}@NE'.format(vxi), bc['CE']+((fc['NE']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE'.format(vxi), bc['HE2']+((fc['HE']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), bc['HE3']-(bc['HE3']/10)*i).execute()
                change(parm, 'charge', ':{}@CZ'.format(vxi), bc['NZ']+((fc['CZ']-bc['NZ'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ1'.format(vxi), bc['HZ1']-(bc['HZ1']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vxi), bc['HZ2']-(bc['HZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ3'.format(vxi), bc['HZ3']-(bc['HZ3']/10)*i).execute()
                change(parm, 'charge', ':{}@NH1'.format(vxi), (fc['NH1']/10)*i).execute()
                change(parm, 'charge', ':{}@HH11'.format(vxi), (fc['HH11']/10)*i).execute()
                change(parm, 'charge', ':{}@HH12'.format(vxi), (fc['HH12']/10)*i).execute()
                change(parm, 'charge', ':{}@NH2'.format(vxi), (fc['NH2']/10)*i).execute()
                change(parm, 'charge', ':{}@HH21'.format(vxi), (fc['HH21']/10)*i).execute()
                change(parm, 'charge', ':{}@HH22'.format(vxi), (fc['HH22']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        NZ = struct.residue_dict[aa].atom_dict['NZ']
        HZ1 = struct.residue_dict[aa].atom_dict['HZ1']
        HZ2 = struct.residue_dict[aa].atom_dict['HZ2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('NE')) 
			elif atom.get_name() == 'HE2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HE')) 
			elif atom.get_name() == 'NZ' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ')) 
			elif atom.get_name() == 'HZ3' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('NH1', NZ, HZ1)) 
				pdb.write(atom.superimposed1('HH11', HZ1)) 
				pdb.write(atom.superimposed2('HH12', HZ1)) 
				pdb.write(atom.halfway_between('NH2', NZ, HZ2)) 
				pdb.write(atom.superimposed1('HH21', HZ2)) 
				pdb.write(atom.superimposed2('HH22', HZ2)) 
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
	he2 = var[2]
	he3 = var[3]
	cz = var[4]
	hz1 = var[5]
	hz2 = var[6]
	hz3 = var[7]
	nh1 = var[8]
	hh1 = var[9]
	nh2 = var[10]
	hh2 = var[11]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ARG-LYS.pdb\n"%vxi)
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
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "C"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "N"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "H"\n'%vxi)
	ctrl.write('set %s.1.24 element "N"\n'%vxi)
	ctrl.write('set %s.1.25 element "H"\n'%vxi)
	ctrl.write('set %s.1.26 element "H"\n'%vxi)
	ctrl.write('set %s.1.27 element "C"\n'%vxi)
	ctrl.write('set %s.1.28 element "O"\n'%vxi)
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
	ctrl.write('set %s.1.16 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.17 name "CZ"\n'%vxi)
	ctrl.write('set %s.1.18 name "HZ1"\n'%vxi)
	ctrl.write('set %s.1.19 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.20 name "HZ3"\n'%vxi)
	ctrl.write('set %s.1.21 name "NH1"\n'%vxi)
	ctrl.write('set %s.1.22 name "HH11"\n'%vxi)
	ctrl.write('set %s.1.23 name "HH12"\n'%vxi)
	ctrl.write('set %s.1.24 name "NH2"\n'%vxi)
	ctrl.write('set %s.1.25 name "HH21"\n'%vxi)
	ctrl.write('set %s.1.26 name "HH22"\n'%vxi)
	ctrl.write('set %s.1.27 name "C"\n'%vxi)
	ctrl.write('set %s.1.28 name "O"\n'%vxi)
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
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, he2))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, he3))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, cz))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, hz1))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, hz2))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, hz3))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, nh1))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, hh1))
	ctrl.write('set %s.1.23 type "%s"\n'%(vxi, hh1))
	ctrl.write('set %s.1.24 type "%s"\n'%(vxi, nh2))
	ctrl.write('set %s.1.25 type "%s"\n'%(vxi, hh2))
	ctrl.write('set %s.1.26 type "%s"\n'%(vxi, hh2))
	ctrl.write('set %s.1.27 type "C"\n'%vxi)
	ctrl.write('set %s.1.28 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.27\n'%(vxi, vxi))
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
	ctrl.write('bond %s.1.14 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.24 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.24 %s.1.26\n'%(vxi, vxi))
	ctrl.write('bond %s.1.27 %s.1.28\n'%(vxi, vxi))
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
	he2 = var[2]
	he3 = var[3]
	cz = var[4]
	hz1 = var[5]
	hz2 = var[6]
	hz3 = var[7]
	nh1 = var[8]
	hh1 = var[9]
	nh2 = var[10]
	hh2 = var[11]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(hd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ne, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(he2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(he3, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cz, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(hz1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hz2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hz3, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd1, lac(p['H1'][0], p['HC'][0], i), lac(p['H1'][1], p['HC'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ne, lac(p['NA'][0], p['CT'][0], i), lac(p['NA'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he2, lac(p['H'][0], p['HP'][0], i), lac(p['H'][1], p['HP'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he3, lac(p['0_H'][0], p['HP'][0], i), lac(p['0_H'][1], p['HP'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz, lac(p['CA'][0], p['NA'][0], i), lac(p['CA'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz1, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz2, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz3, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), nh1, lac(p['NA'][0], p['0_N'][0], i), lac(p['NA'][1], p['0_N'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh1, lac(p['H'][0], p['0_H'][0], i), lac(p['H'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), nh2, lac(p['NA'][0], p['0_N'][0], i), lac(p['NA'][1], p['0_N'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh2, lac(p['H'][0], p['0_H'][0], i), lac(p['H'][1], p['0_H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hd1), lac(p['CT_HC'][0], p['CT_HC'][0], i), lac(p['CT_HC'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', ne), lac(p['CT_N'][0], p['CT_CT'][0], i), lac(p['CT_N'][1], p['CT_CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne, he2), lac(p['NA_H'][0], p['CT_HC'][0], i), lac(p['NA_H'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne, he3), lac(p['CT_sC'][0], p['CT_HC'][0], i), lac(p['CT_sC'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne, cz), lac(p['CA_N'][0], p['C_N3'][0], i), lac(p['CA_N'][1], p['C_N3'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, hz1), lac(p['NA_sN'][0], p['NA_H'][0], i), lac(p['NA_sN'][1], p['NA_sN'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, hz2), lac(p['NA_sN'][0], p['NA_H'][0], i), lac(p['NA_sN'][1], p['NA_sN'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, hz3), lac(p['NA_sC'][0], p['NA_H'][0], i), lac(p['NA_sC'][1], p['NA_sN'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, nh1), lac(p['CA_N'][0], p['NA_mH2'][0], i), lac(p['CA_N'][1], p['NA_mH2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, nh2), lac(p['CA_N'][0], p['NA_mH2'][0], i), lac(p['CA_N'][1], p['NA_mH2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(nh1, hh1), lac(p['NA_H'][0], p['NA_mH'][0], i), lac(p['NA_H'][1], p['NA_mH'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(nh2, hh2), lac(p['NA_H'][0], p['NA_mH'][0], i), lac(p['NA_H'][1], p['NA_mH'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hd1), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, 'CT', hd1), lac(p['H_C_H'][0], p['H_C_H'][0], i), lac(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', ne), lac(p['C_C_N2'][0], p['C_C_C'][0], i), lac(p['C_C_N2'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, 'CT', ne), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ne, he2), lac(p['C_N_H'][0], p['C_C_H'][0], i), lac(p['C_N_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ne, he3), lac(p['C_N_C'][0], p['C_C_H'][0], i), lac(p['C_N_C'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', ne, cz), lac(p['C_N_C'][0], p['C_C_N'][0], i), lac(p['C_C_N'][1], p['C_C_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he2, ne, cz), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he3, ne, cz), lac(p['Close'][0], p['C_C_H'][0], i), lac(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he2, ne, he3), lac(p['F_CA_CA_HA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, hz1), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, hz2), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, hz3), lac(p['Close'][0], p['C_C_H'][0], i), lac(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, nh1), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne, cz, nh2), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz1, cz, nh1), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz1, cz, nh2), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz2, cz, nh1), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz2, cz, nh2), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz3, cz, nh1), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz3, cz, nh2), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(nh1, cz, nh2), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz, nh1, hh1), lac(p['F_CA_CA_CA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_CA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz, nh2, hh2), lac(p['F_CA_CA_CA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_CA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh1, nh1, hh1), lac(p['H_N_H'][0], p['Close'][0], i), lac(p['H_N_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh2, nh2, hh2), lac(p['H_N_H'][0], p['Close'][0], i), lac(p['H_N_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz1, cz, hz2), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz2, cz, hz3), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz1, cz, hz3), lac(p['F_CA_CA_CA'][0], p['H_C_H'][0], i), lac(p['F_CA_CA_CA'][1], p['H_C_H'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ne, he2), lac(p['X_CT_N_X'][0], p['C_C_C_H'][0], i), lac(p['X_CT_N_X'][1], p['C_C_C_H'][1], i), lac(p['X_CT_N_X'][2], p['C_C_C_H'][2], i), lac(p['X_CT_N_X'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ne, he3), lac(p['0_1'][0], p['C_C_C_H'][0], i), lac(p['0_1'][1], p['C_C_C_H'][1], i), lac(p['0_1'][2], p['C_C_C_H'][2], i), lac(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, 'CT', ne, he2), lac(p['X_CT_N_X'][0], p['H_C_C_H'][0], i), lac(p['X_CT_N_X'][1], p['H_C_C_H'][1], i), lac(p['X_CT_N_X'][2], p['H_C_C_H'][2], i), lac(p['X_CT_N_X'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, 'CT', ne, he3), lac(p['0_1'][0], p['H_C_C_H'][0], i), lac(p['0_1'][1], p['H_C_C_H'][1], i), lac(p['0_1'][2], p['H_C_C_H'][2], i), lac(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', ne, cz), lac(p['X_CT_N_X'][0], p['X_C_C_X'][0], i), lac(p['X_CT_N_X'][1], p['X_C_C_X'][1], i), lac(p['X_CT_N_X'][2], p['X_C_C_X'][2], i), lac(p['X_CT_N_X'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, 'CT', ne, cz), lac(p['X_CT_N_X'][0], p['X_C_C_X'][0], i), lac(p['X_CT_N_X'][1], p['X_C_C_X'][1], i), lac(p['X_CT_N_X'][2], p['X_C_C_X'][2], i), lac(p['X_CT_N_X'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, hz1), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, hz2), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, hz3), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he2, ne, cz, hz1), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he2, ne, cz, hz2), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he2, ne, cz, hz3), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, ne, cz, hz1), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, ne, cz, hz2), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, ne, cz, hz3), lac(p['0_1'][0], p['X_C_C_X'][0], i), lac(p['0_1'][1], p['X_C_C_X'][1], i), lac(p['0_1'][2], p['X_C_C_X'][2], i), lac(p['0_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, nh1), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, nh2), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he2, ne, cz, nh1), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he2, ne, cz, nh2), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, ne, cz, nh1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, ne, cz, nh2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne, cz, nh1, hh1), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne, cz, nh2, hh2), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz1, cz, nh1, hh1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz1, cz, nh2, hh2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz2, cz, nh1, hh1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz2, cz, nh2, hh2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz3, cz, nh1, hh1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz3, cz, nh2, hh2), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nh2, cz, nh1, hh1), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(nh1, cz, nh2, hh2), lac(p['X_CA_N_X'][0], p['Ring_0'][0], i), lac(p['X_CA_N_X'][1], p['Ring_0'][1], i), lac(p['X_CA_N_X'][2], p['Ring_0'][2], i), lac(p['X_CA_N_X'][3], p['Ring_0'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz, hh1, nh1, hh1), lac(p['Ami_imp'][0], p['Imp_0'][0], i), lac(p['Ami_imp'][1], p['Imp_0'][1], i), lac(p['Ami_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz, hh2, nh2, hh2), lac(p['Ami_imp'][0], p['Imp_0'][0], i), lac(p['Ami_imp'][1], p['Imp_0'][1], i), lac(p['Ami_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne, nh1, cz, nh2), lac(p['Car_imp'][0], p['Imp_0'][0], i), lac(p['Car_imp'][1], p['Imp_0'][1], i), lac(p['Car_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', ne, cz, he2), lac(p['Car_imp'][0], p['Imp_0'][0], i), lac(p['Car_imp'][1], p['Imp_0'][1], i), lac(p['Car_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd1, lac(p['H1'][2], p['HC'][2], i), lac(p['H1'][3], p['HC'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ne, lac(p['NA'][2], p['CT'][2], i), lac(p['NA'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he2, lac(p['H'][2], p['HP'][2], i), lac(p['H'][3], p['HP'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he3, lac(p['0_H'][2], p['HP'][2], i), lac(p['0_H'][3], p['HP'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz, lac(p['CA'][2], p['NA'][2], i), lac(p['CA'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz1, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz2, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz3, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), nh1, lac(p['NA'][2], p['0_N'][2], i), lac(p['NA'][3], p['0_N'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh1, lac(p['H'][2], p['0_H'][2], i), lac(p['H'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), nh2, lac(p['NA'][2], p['0_N'][2], i), lac(p['NA'][3], p['0_N'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh2, lac(p['H'][2], p['0_H'][2], i), lac(p['H'][3], p['0_H'][3], i))
