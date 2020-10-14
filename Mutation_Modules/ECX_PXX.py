# ECX to PXX Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/ECX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/PXX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
#		changeLJPair(parm, ':{}@HD2'.format(vxi), ':{}@HZ1'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HD2'.format(vxi), ':{}@HD'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HZ1'.format(vxi), ':{}@HD'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HZ1'.format(vxi), ':{}@CK'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@CZ1'.format(vxi), ':{}@CK'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@CZ1'.format(vxi), ':{}@CH2'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@CZ2'.format(vxi), ':{}@CH1'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HZ1'.format(vxi), ':{}@HH1'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HZ1'.format(vxi), ':{}@HH2'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HZ1'.format(vxi), ':{}@HZ2'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HH1'.format(vxi), ':{}@HH2'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HH1'.format(vxi), ':{}@HZ2'.format(vxi), '0', '0').execute()
#		changeLJPair(parm, ':{}@HH2'.format(vxi), ':{}@HZ2'.format(vxi), '0', '0').execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@SG'.format(vxi), bc['SG']+((fc['SG']-bc['SG'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), (fc['CD']/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), (fc['HD2']/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), (fc['HD3']/10)*i).execute()
                change(parm, 'charge', ':{}@CE'.format(vxi), (fc['CE']/10)*i).execute()
                change(parm, 'charge', ':{}@CZ1'.format(vxi), bc['CD']+((fc['CZ1']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ1'.format(vxi), bc['HD2']+((fc['HZ1']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD'.format(vxi), bc['HD3']-(bc['HD3']/10)*i).execute()
                change(parm, 'charge', ':{}@CH1'.format(vxi), (fc['CH1']/10)*i).execute()
                change(parm, 'charge', ':{}@HH1'.format(vxi), (fc['HH1']/10)*i).execute()
                change(parm, 'charge', ':{}@CT'.format(vxi), (fc['CT']/10)*i).execute()
                change(parm, 'charge', ':{}@CH2'.format(vxi), (fc['CH2']/10)*i).execute()
                change(parm, 'charge', ':{}@HH2'.format(vxi), (fc['HH2']/10)*i).execute()
                change(parm, 'charge', ':{}@CZ2'.format(vxi), (fc['CZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vxi), (fc['HZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@CK'.format(vxi), bc['CE']+((fc['CK']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HK2'.format(vxi), bc['HE2']+((fc['HK2']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HK3'.format(vxi), bc['HE3']+((fc['HK3']-bc['HE3'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        SG = struct.residue_dict[aa].atom_dict['SG']
        CD = struct.residue_dict[aa].atom_dict['CD']
        CE = struct.residue_dict[aa].atom_dict['CE']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'SG' and res.get_resname() == vxi:
				pdb.write(atom.formatted())
				pdb.write(atom.thirdone_between('CD', SG, CD)) 
				pdb.write(atom.superimposed1('HD2', SG)) 
				pdb.write(atom.superimposed2('HD3', SG)) 
				pdb.write(atom.thirdtwo_between('CE', SG, CD)) 
			elif atom.get_name() == 'CD' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ1')) 
			elif atom.get_name() == 'HD2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HZ1')) 
			elif atom.get_name() == 'HD3' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HD')) 
				pdb.write(atom.thirdone_between1('CH1', CD, CE)) 
				pdb.write(atom.superimposed1('HH1', CD)) 
				pdb.write(atom.thirdtwo_between('CT', CD, CE)) 
				pdb.write(atom.thirdone_between1('CH2', CD, CE)) 
				pdb.write(atom.thirdtwo_between1('HH2', CD, CE)) 
				pdb.write(atom.superimposed2('CZ2', CD)) 
				pdb.write(atom.thirdtwo_between1('HZ2', SG, CD)) 
			elif atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CK')) 
			elif atom.get_name() == 'HE2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HK2')) 
			elif atom.get_name() == 'HE3' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HK3')) 
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
	cd = var[0]
	hd1 = var[1]
	ce = var[2]
	cz1 = var[3]
	hz1 = var[4]
	hd = var[5]
	ch1 = var[6]
	hh1 = var[7]
	ct = var[8]
	ch2 = var[9]
	hh2 = var[10]
	cz2 = var[11]
	hz2 = var[12]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ECX-PXX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "S"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "C"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "C"\n'%vxi)
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
	ctrl.write('set %s.1.8 name "SG"\n'%vxi)
	ctrl.write('set %s.1.9 name "CD"\n'%vxi)
	ctrl.write('set %s.1.10 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.12 name "CE"\n'%vxi)
	ctrl.write('set %s.1.13 name "CZ1"\n'%vxi)
	ctrl.write('set %s.1.14 name "HZ1"\n'%vxi)
	ctrl.write('set %s.1.15 name "HD"\n'%vxi)
	ctrl.write('set %s.1.16 name "CH1"\n'%vxi)
	ctrl.write('set %s.1.17 name "HH1"\n'%vxi)
	ctrl.write('set %s.1.18 name "CT"\n'%vxi)
	ctrl.write('set %s.1.19 name "CH2"\n'%vxi)
	ctrl.write('set %s.1.20 name "HH2"\n'%vxi)
	ctrl.write('set %s.1.21 name "CZ2"\n'%vxi)
	ctrl.write('set %s.1.22 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.23 name "CK"\n'%vxi)
	ctrl.write('set %s.1.24 name "HK2"\n'%vxi)
	ctrl.write('set %s.1.25 name "HK3"\n'%vxi)
	ctrl.write('set %s.1.26 name "C"\n'%vxi)
	ctrl.write('set %s.1.27 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "H1"\n'%vxi)
	ctrl.write('set %s.1.7 type "H1"\n'%vxi)
	ctrl.write('set %s.1.8 type "S"\n'%vxi)
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, cd))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, hd1))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, hd1))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, ce))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, cz1))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, hz1))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, hd))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, ch1))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, hh1))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, ct))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, ch2))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, hh2))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, cz2))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, hz2))
	ctrl.write('set %s.1.23 type "CT"\n'%vxi)
	ctrl.write('set %s.1.24 type "H1"\n'%vxi)
	ctrl.write('set %s.1.25 type "H1"\n'%vxi)
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
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.26 %s.1.27\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect2 %s.1.CK\n'%(vxi, vxi))
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
	cd = var[0]
	hd1 = var[1]
	ce = var[2]
	cz1 = var[3]
	hz1 = var[4]
	hd = var[5]
	ch1 = var[6]
	hh1 = var[7]
	ct = var[8]
	ch2 = var[9]
	hh2 = var[10]
	cz2 = var[11]
	hz2 = var[12]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(cd, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ce, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(cz1, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hz1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ch1, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hh1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ct, 'C', 'sp2')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd, cal(p['0_C'][0], p['CT'][0], i), cal(p['0_C'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd1, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ce, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz1, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz1, cal(p['H1'][0], p['HA'][0], i), cal(p['H1'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd, cal(p['H1'][0], p['0_H'][0], i), cal(p['H1'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ch1, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh1, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ct, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ch2, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh2, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz2, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz2, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('S ', cd), cal(p['S_tCT'][0], p['CT_S'][0], i), cal(p['S_tCT'][1], p['CT_S'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd, hd1), cal(p['CT_tHC'][0], p['CT_HC'][0], i), cal(p['CT_tHC'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd, ce), cal(p['CT_tCA'][0], p['CT_CA'][0], i), cal(p['CT_tCA'][1], p['CT_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce, cz1), cal(p['CT_tCA'][0], p['CA_CA'][0], i), cal(p['CT_tCA'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce, cz2), cal(p['CT_tCA'][0], p['CA_CA'][0], i), cal(p['CT_tCA'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz1, hz1), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz1, hd), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz1, ch1), cal(p['CT_tCT'][0], p['CA_CA'][0], i), cal(p['CT_tCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch1, hh1), cal(p['CT_tHA'][0], p['CA_HA'][0], i), cal(p['CT_tHA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch1, ct), cal(p['CT_tCT'][0], p['CA_CA'][0], i), cal(p['CT_tCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ct, 'CT'), cal(p['CT_tCT'][0], p['CT_CA'][0], i), cal(p['CT_tCT'][1], p['CT_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ct, ch2), cal(p['CT_tCT'][0], p['CA_CA'][0], i), cal(p['CT_tCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch2, hh2), cal(p['CT_tHA'][0], p['CA_HA'][0], i), cal(p['CT_tHA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch2, cz2), cal(p['CT_tCT'][0], p['CA_CA'][0], i), cal(p['CT_tCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz2, hz2), cal(p['CT_tHA'][0], p['CA_HA'][0], i), cal(p['CT_tHA'][1], p['CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'S ', cd), cal(p['C_S_C'][0], p['C_S_C'][0], i), cal(p['C_S_C'][1], p['C_S_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', cd, hd1), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, cd, hd1), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', cd, ce), cal(p['Dritt'][0], p['CA_C_S'][0], i), cal(p['Dritt'][1], p['CA_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd1, cd, ce), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd, ce, cz1), cal(p['Dritt'][0], p['F_CT_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd, ce, cz2), cal(p['Dritt'][0], p['F_CT_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz1, ce, cz2), cal(p['Close'][0], p['F_CT_CA_CA'][0], i), cal(p['Close'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce, cz1, hz1), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce, cz1, hd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce, cz1, ch1), cal(p['C_C_S'][0], p['F_CA_CA_CA'][0], i), cal(p['C_C_S'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce, cz2, hz2), cal(p['C_C_S'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_S'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz1, cz1, hd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz1, cz1, ch1), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd, cz1, ch1), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz1, ch1, hh1), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz1, ch1, ct), cal(p['Dritt'][0], p['F_CA_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh1, ch1, ct), cal(p['Dritt'][0], p['F_CA_CA_HA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch1, ct, 'CT'), cal(p['Dritt'][0], p['F_CT_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch2, ct, 'CT'), cal(p['Dritt'][0], p['F_CT_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ct, 'CT', 'H1'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ct, 'CT', 'S '), cal(p['C_C_S'][0], p['CA_C_S'][0], i), cal(p['C_C_S'][1], p['CA_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch1, ct, ch2), cal(p['Close'][0], p['F_CA_CA_CA'][0], i), cal(p['Close'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ct, ch2, hh2), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ct, ch2, cz2), cal(p['Dritt'][0], p['F_CA_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hh2, ch2, cz2), cal(p['Dritt'][0], p['F_CA_CA_HA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch2, cz2, hz2), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch2, cz2, ce), cal(p['C_C_S'][0], p['F_CA_CA_CA'][0], i), cal(p['C_C_S'][1], p['F_CA_CA_CA'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'S ', cd, hd1), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'S ', cd, ce), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', cd, ce, cz1), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', cd, ce, cz2), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, cd, ce, cz1), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, cd, ce, cz2), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, ce, cz1, ch1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, ce, cz1, hz1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, ce, cz1, hd), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, ce, cz2, hz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd, ce, cz2, ch2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce, cz1, ch1, hh1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz1, cz1, ch1, hh1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd, cz1, ch1, hh1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce, cz1, ch1, ct), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz1, cz1, ch1, ct), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd, cz1, ch1, ct), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz1, ch1, ct, ch2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch1, ct, ch2, cz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch1, ct, ch2, hh2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ct, ch2, cz2, ce), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh2, ch2, cz2, ce), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ct, ch2, cz2, hz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh2, ch2, cz2, hz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh1, ch1, ct, ch2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz2, cz2, ce, cz1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch2, cz2, ce, cz1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz1, cz1, ce, cz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd, cz1, ce, cz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch1, cz1, ce, cz2), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz1, ch1, ct, 'CT'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh1, ch1, ct, 'CT'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz2, ch2, ct, 'CT'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh2, ch2, ct, 'CT'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch1, ct, 'CT', 'H1'), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch1, ct, 'CT', 'S '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch2, ct, 'CT', 'H1'), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ch2, ct, 'CT', 'S '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		#Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X', 'X', cd, hd1), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd, cal(p['0_C'][2], p['0_C'][2], i), cal(p['0_C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd1, cal(p['0_H'][2], p['0_H'][2], i), cal(p['0_H'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ce, cal(p['0_C'][2], p['0_C'][2], i), cal(p['0_C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz1, cal(p['0_C'][2], p['0_C'][2], i), cal(p['0_C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz1, cal(p['0_H'][2], p['0_H'][2], i), cal(p['0_H'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd, cal(p['0_H'][2], p['0_H'][2], i), cal(p['0_H'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ch1, cal(p['0_C'][2], p['0_C'][2], i), cal(p['0_C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh1, cal(p['0_H'][2], p['0_H'][2], i), cal(p['0_H'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ct, cal(p['0_C'][2], p['0_C'][2], i), cal(p['0_C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ch2, cal(p['0_C'][2], p['0_C'][2], i), cal(p['0_C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh2, cal(p['0_H'][2], p['0_H'][2], i), cal(p['0_H'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz2, cal(p['0_C'][2], p['0_C'][2], i), cal(p['0_C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz2, cal(p['0_H'][2], p['0_H'][2], i), cal(p['0_H'][3], p['0_H'][3], i))

#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd1, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ce, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz1, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz1, cal(p['H1'][2], p['HA'][2], i), cal(p['H1'][3], p['HA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd, cal(p['H1'][2], p['0_H'][2], i), cal(p['H1'][3], p['0_H'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ch1, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh1, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ct, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ch2, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh2, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz2, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
#		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz2, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
