# TYR to TRP Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
	bc = {}
        with open('Param_files/AminoAcid/TYR.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/TRP.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                #changeLJPair(parm, ':{}@HH2 :{}@HZ1 0 0'.format(vxi, vxi)).execute()
                #changeLJPair(parm, ':{}@HE3 :{}@HZ2 0 0'.format(vxi, vxi)).execute()
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
                change(parm, 'charge', ':{}@NE1'.format(vxi), bc['CE1']+((fc['NE1']-bc['CE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), bc['HE1']+((fc['HE1']-bc['HE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE3'.format(vxi), bc['HD2']+((fc['CE3']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), (fc['HE3']/10)*i).execute()
                change(parm, 'charge', ':{}@CZ'.format(vxi), bc['CZ']-(bc['CZ']/10)*i).execute()
                change(parm, 'charge', ':{}@OH'.format(vxi), bc['OH']-(bc['OH']/10)*i).execute()
                change(parm, 'charge', ':{}@HH'.format(vxi), bc['HH']-(bc['HH']/10)*i).execute()
                change(parm, 'charge', ':{}@CZ2'.format(vxi), bc['HE2']+((fc['CZ1']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ3'.format(vxi), (fc['CZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ3'.format(vxi), (fc['HZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vxi), (fc['HZ1']/10)*i).execute()
                change(parm, 'charge', ':{}@CH2'.format(vxi), (fc['CH']/10)*i).execute()
                change(parm, 'charge', ':{}@HH2'.format(vxi), (fc['HH']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	CD2 = struct.residue_dict[aa].atom_dict['CD2']
	HD2 = struct.residue_dict[aa].atom_dict['HD2']
	CE2 = struct.residue_dict[aa].atom_dict['CE2']
	HE2 = struct.residue_dict[aa].atom_dict['HE2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'CE1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('NE1'))
			elif atom.get_name() == 'HD2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CE3'))
				pdb.write(atom.superimposed1('HE3', CD2))
			elif atom.get_name() == 'HE2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ2'))
				pdb.write(atom.superimposed1('HZ2', CE2))
				pdb.write(atom.thirdone_between('CZ3', HD2, HE2))
				pdb.write(atom.superimposed1('HZ3', HD2))
				pdb.write(atom.thirdtwo_between('CH2', HD2, HE2))
				pdb.write(atom.superimposed1('HH2', HE2))
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

def makevxi_alt(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	CD1 = struct.residue_dict[aa].atom_dict['CD1']
	CD2 = struct.residue_dict[aa].atom_dict['CD2']
	HD1 = struct.residue_dict[aa].atom_dict['HD1']
	CE = struct.residue_dict[aa].atom_dict['CE']
	HE = struct.residue_dict[aa].atom_dict['HE']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'CD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CD2'))
			elif atom.get_name() == 'HD1' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CE3'))
			elif atom.get_name() == 'CD2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CD1'))
				pdb.write(atom.halfway_between('NE1', CD2, CE))
				pdb.write(atom.superimposed1('HE1', CD2))
			elif atom.get_name() == 'HD2' and res.get_resname() == vxi:
				pdb.write(atom.change_name('HD1'))
				pdb.write(atom.superimposed1('HE3', CD1))
			elif atom.get_name() == 'CE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CE2'))
			elif atom.get_name() == 'HE' and res.get_resname() == vxi:
				pdb.write(atom.change_name('CZ2'))
				pdb.write(atom.superimposed1('HZ2', CE))
				pdb.write(atom.thirdone_between('CZ3', HD1, HE))
				pdb.write(atom.superimposed1('HZ3', HD1))
				pdb.write(atom.thirdtwo_between('CH2', HD1, HE))
				pdb.write(atom.superimposed1('HH2', HE))
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

def lib_make(ff, outputfile, vxi='VXI', cg='1c', cd1='2c', hd1='2h', ne1='3n', he1='3h', cz='0c', oh='0o', hh='0h', cd2='4c', ce2='5c', ce3='6c', he3='6h', cz1='7c', hz1='7h', cz2='8c', hz2='8h', ch='9c', hh2='9h'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/TYR-TRP.pdb\n"%vxi)
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
	ctrl.write('set %s.1.12 element "N"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "O"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "C"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "C"\n'%vxi)
	ctrl.write('set %s.1.23 element "H"\n'%vxi)
	ctrl.write('set %s.1.24 element "C"\n'%vxi)
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
	ctrl.write('set %s.1.9 name "CD1"\n'%vxi)
	ctrl.write('set %s.1.10 name "HD1"\n'%vxi)
	ctrl.write('set %s.1.11 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.12 name "NE1"\n'%vxi)
	ctrl.write('set %s.1.13 name "HE1"\n'%vxi)
	ctrl.write('set %s.1.14 name "CE3"\n'%vxi)
	ctrl.write('set %s.1.15 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.16 name "CZ"\n'%vxi)
	ctrl.write('set %s.1.17 name "OH"\n'%vxi)
	ctrl.write('set %s.1.18 name "HH"\n'%vxi)
	ctrl.write('set %s.1.19 name "CE2"\n'%vxi)
	ctrl.write('set %s.1.20 name "CZ2"\n'%vxi)
	ctrl.write('set %s.1.21 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.22 name "CZ3"\n'%vxi)
	ctrl.write('set %s.1.23 name "HZ3"\n'%vxi)
	ctrl.write('set %s.1.24 name "CH2"\n'%vxi)
	ctrl.write('set %s.1.25 name "HH2"\n'%vxi)
	ctrl.write('set %s.1.26 name "C"\n'%vxi)
	ctrl.write('set %s.1.27 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, cg))
	ctrl.write('set %s.1.9  type "%s"\n'%(vxi, cd1))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, hd1))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, cd2))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, ne1))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, he1))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, ce3))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, he3))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, cz))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, oh))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, hh))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, ce2))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, cz1))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, hz1))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, cz2))
	ctrl.write('set %s.1.23 type "%s"\n'%(vxi, hz2))
	ctrl.write('set %s.1.24 type "%s"\n'%(vxi, ch))
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
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.22 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.22 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.24 %s.1.25\n'%(vxi, vxi))
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

def stock_add_to_all(cg='1c', cd1='2c', hd1='2h', ne1='3n', he1='3h', cd2='4c', ce2='5c', cz='0c', oh='0o', hh='0h', ce3='6c', he3='6h', cz1='7c', hz1='7h', cz2='8c', hz2='8h', ch='9c', hh2='9h'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(cg, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(cd1, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ne1, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(he1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cd2, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(ce2, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(ce3, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(he3, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cz, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(oh, 'O', 'sp3')
	Frcmod_creator.TYPE_insert(hh, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cz1, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hz1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(cz2, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hz2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(ch, 'C', 'sp2')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cg, lac(p['CA'][0], p['CA'][0], i), lac(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd1, lac(p['CA'][0], p['CA'][0], i), lac(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd1, lac(p['H4'][0], p['HA'][0], i), lac(p['H4'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ne1, lac(p['NA'][0], p['CA'][0], i), lac(p['NA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he1, lac(p['H'][0], p['HA'][0], i), lac(p['H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ce2, lac(p['CA'][0], p['CA'][0], i), lac(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd2, lac(p['CA'][0], p['CA'][0], i), lac(p['CA'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ce3, lac(p['CA'][0], p['HA'][0], i), lac(p['CA'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), he3, lac(p['HA'][0], p['0_H'][0], i), lac(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz, lac(p['0_C'][0], p['CA'][0], i), lac(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), oh, lac(p['0_O'][0], p['OH'][0], i), lac(p['0_O'][1], p['OH'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh, lac(p['0_H'][0], p['HO'][0], i), lac(p['0_H'][1], p['HO'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz1, lac(p['CA'][0], p['HA'][0], i), lac(p['CA'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz1, lac(p['HA'][0], p['0_H'][0], i), lac(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cz2, lac(p['CA'][0], p['0_C'][0], i), lac(p['CA'][1], p['0_C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hz2, lac(p['HA'][0], p['0_H'][0], i), lac(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ch, lac(p['CA'][0], p['0_C'][0], i), lac(p['CA'][1], p['0_C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hh2, lac(p['HA'][0], p['0_H'][0], i), lac(p['HA'][1], p['0_H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', cg), lac(p['CT_C*'][0], p['CT_CA'][0], i), lac(p['CT_C*'][1], p['CT_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cg, cd1), lac(p['C*_CW'][0], p['CA_CA'][0], i), lac(p['C*_CW'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cg, cd2), lac(p['C*_CB'][0], p['CA_CA'][0], i), lac(p['C*_CB'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce2, cd2), lac(p['CB_CN'][0], p['CA_CA'][0], i), lac(p['CB_CN'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd1, hd1), lac(p['CA_HA'][0], p['CA_HA'][0], i), lac(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd1, ne1), lac(p['CW_NA'][0], p['CA_CA'][0], i), lac(p['CW_NA'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne1, he1), lac(p['NA_H'][0], p['CA_HA'][0], i), lac(p['NA_H'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ne1, cz), lac(p['CA_mCN'][0], p['CA_CA'][0], i), lac(p['CA_mCN'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, oh), lac(p['OH_mCN'][0], p['CA_OH'][0], i), lac(p['OH_mCN'][1], p['CA_OH'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(oh, hh), lac(p['HO_mNA'][0], p['OH_HO'][0], i), lac(p['HO_mNA'][1], p['OH_HO'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz, ce2), lac(p['CA_mCN'][0], p['CA_CA'][0], i), lac(p['CA_mCN'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce2, cz1), lac(p['CA_CA'][0], p['CA_HA'][0], i), lac(p['CA_CA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd2, ce3), lac(p['CA_CB'][0], p['CA_HA'][0], i), lac(p['CA_CB'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce3, he3), lac(p['CA_HA'][0], p['CA_HA'][0], i), lac(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz1, hz1), lac(p['CA_HA'][0], p['CA_HA'][0], i), lac(p['CA_HA'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce3, cz2), lac(p['CA_CA'][0], p['CA_tCA2'][0], i), lac(p['CA_CA'][1], p['CA_tCA2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz2, ch), lac(p['CA_CA'][0], p['CA_tCA2'][0], i), lac(p['CA_CA'][1], p['CA_tCA2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz1, ch), lac(p['CA_CA'][0], p['CA_tCA2'][0], i), lac(p['CA_CA'][1], p['CA_tCA2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cz2, hz2), lac(p['CA_HA'][0], p['HA_tCA2'][0], i), lac(p['CA_HA'][1], p['HA_tCA2'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ch, hh2), lac(p['CA_HA'][0], p['HA_tCA2'][0], i), lac(p['CA_HA'][1], p['HA_tCA2'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', cg), lac(p['C_C_C*'][0], p['C_C_CA'][0], i), lac(p['C_C_C*'][1], p['C_C_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', cg), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', cg, cd1), lac(p['CT_C*_CW'][0], p['F_CT_CA_CA'][0], i), lac(p['CT_C*_CW'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', cg, cd2), lac(p['CT_C*_CB'][0], p['F_CT_CA_CA'][0], i), lac(p['CT_C*_CB'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd1, cg, cd2), lac(p['CW_C*_CB'][0], p['F_CA_CA_CA'][0], i), lac(p['CW_C*_CB'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cg, cd1, ne1), lac(p['C*_CW_NA'][0], p['F_CA_CA_CA'][0], i), lac(p['C*_CW_NA'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cg, cd1, hd1), lac(p['F_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd1, ne1, cz), lac(p['CW_NA_CN'][0], p['F_CA_CA_CA'][0], i), lac(p['CW_NA_CN'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he1, ne1, cz), lac(p['CW_NA_H'][0], p['F_CA_CA_HA'][0], i), lac(p['CW_NA_CN'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne1, cz, ce2), lac(p['Dritt'][0], p['F_CA_CA_CA'][0], i), lac(p['Dritt'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne1, cz, oh), lac(p['Close'][0], p['F_CA_CA_HA'][0], i), lac(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce2, cz, oh), lac(p['Dritt'][0], p['F_CA_CA_HA'][0], i), lac(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz, oh, hh), lac(p['CW_NA_H'][0], p['CA_O_H'][0], i), lac(p['CW_NA_H'][1], p['CA_O_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ne1, cd1, hd1), lac(p['F_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), lac(p['F_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(he1, ne1, cd1), cal(p['F_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['F_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cg, cd2, ce3), lac(p['C*_CB_CA'][0], p['F_CA_CA_HA'][0], i), lac(p['C*_CB_CA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cg, cd2, ce2), lac(p['C*_CB_CN'][0], p['F_CA_CA_CA'][0], i), lac(p['C*_CB_CN'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd2, ce2, cz), lac(p['CB_CN_NA'][0], p['F_CA_CA_CA'][0], i), lac(p['CB_CN_NA'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd2, ce3, he3), lac(p['F_CA_CA_HA'][0], p['Close'][0], i), lac(p['F_CA_CA_HA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce2, cd2, ce3), lac(p['CA_CB_CN'][0], p['F_CA_CA_HA'][0], i), lac(p['CA_CB_CN'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd2, ce3, cz2), lac(p['F_CA_CA_CA'][0], p['A_Trp'][0], i), lac(p['F_CA_CA_CA'][1], p['A_Trp'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz2, ce3, he3), lac(p['F_CA_CA_HA'][0], p['A_Trp'][0], i), lac(p['F_CA_CA_HA'][1], p['A_Trp'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce2, cz1, hz1), lac(p['F_CA_CA_HA'][0], p['A_Trp'][0], i), lac(p['F_CA_CA_HA'][1], p['A_Trp'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce2, cz1, ch), lac(p['F_CA_CA_CA'][0], p['A_Trp'][0], i), lac(p['F_CA_CA_CA'][1], p['A_Trp'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd2, ce2, cz1), lac(p['CA_CN_CB'][0], p['F_CA_CA_HA'][0], i), lac(p['CA_CN_CB'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz1, ce2, cz), lac(p['CA_CN_NA'][0], p['F_CA_CA_HA'][0], i), lac(p['CA_CN_NA'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz2, cz2, ce3), lac(p['F_CA_CA_HA'][0], p['Close'][0], i), lac(p['F_CA_CA_HA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ch, cz2, ce3), lac(p['F_CA_CA_CA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_CA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz1, ch, cz2), lac(p['F_CA_CA_CA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_CA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz1, cz1, ch), lac(p['F_CA_CA_HA'][0], p['Close'][0], i), lac(p['F_CA_CA_HA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz1, ch, hh2), lac(p['F_CA_CA_HA'][0], p['Close'][0], i), lac(p['F_CA_CA_HA'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cz2, ch, hh2), lac(p['F_CA_CA_HA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_HA'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hz2, cz2, ch), lac(p['F_CA_CA_HA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_HA'][1], p['Dritt'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, cd1), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', cg, cd2), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd1), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', cg, cd2), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd1, ne1), lac(p['X_C_CW_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CW_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CW_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CW_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd1, hd1), lac(p['X_C_CW_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CW_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CW_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CW_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, cg, cd1, ne1), lac(p['X_C_CW_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CW_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CW_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CW_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, cg, cd1, hd1), lac(p['X_C_CW_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CW_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CW_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CW_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd1, ne1, he1), lac(p['X_CW_NA_X'][0], p['Ring_Dihe'][0], i), lac(p['X_CW_NA_X'][1], p['Ring_Dihe'][1], i), lac(p['X_CW_NA_X'][2], p['Ring_Dihe'][2], i), lac(p['X_CW_NA_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd1, ne1, cz), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, cd1, ne1, he1), lac(p['X_CW_NA_X'][0], p['Ring_Dihe'][0], i), lac(p['X_CW_NA_X'][1], p['Ring_Dihe'][1], i), lac(p['X_CW_NA_X'][2], p['Ring_Dihe'][2], i), lac(p['X_CW_NA_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hd1, cd1, ne1, cz), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce2, cz, ne1, he1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(oh, cz, ne1, he1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce2, cz, ne1, cd1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(oh, cz, ne1, cd1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, ce2), lac(p['X_C_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', cg, cd2, ce3), lac(p['X_C_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd1, cg, cd2, ce2), lac(p['X_C_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd1, cg, cd2, ce3), lac(p['X_C_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_C_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_C_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_C_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, ce3, he3), lac(p['X_CA_CB_X'][0], p['Ring_0'][0], i), lac(p['X_CA_CB_X'][1], p['Ring_0'][1], i), lac(p['X_CA_CB_X'][2], p['Ring_0'][2], i), lac(p['X_CA_CB_X'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, ce3, cz2), lac(p['X_CA_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_CA_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_CA_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_CA_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce2, cd2, ce3, he3), lac(p['X_CA_CB_X'][0], p['Ring_0'][0], i), lac(p['X_CA_CB_X'][1], p['Ring_0'][1], i), lac(p['X_CA_CB_X'][2], p['Ring_0'][2], i), lac(p['X_CA_CB_X'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce2, cd2, ce3, cz2), lac(p['X_CA_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_CA_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_CA_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_CA_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, ce2, cz), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cg, cd2, ce2, cz1), lac(p['X_CN_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_CN_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_CN_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_CN_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, ce2, cz, oh), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, ce2, cz, ne1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce3, cd2, ce2, cz), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce3, cd2, ce2, cz1), lac(p['X_CN_CB_X'][0], p['Ring_Dihe'][0], i), lac(p['X_CN_CB_X'][1], p['Ring_Dihe'][1], i), lac(p['X_CN_CB_X'][2], p['Ring_Dihe'][2], i), lac(p['X_CN_CB_X'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, ce2, cz1, hz1), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, ce2, cz1, ch), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh, oh, cz, ce2), lac(p['Ring_0'][0], p['Ring_Dihe'][0], i), lac(p['Ring_0'][1], p['Ring_Dihe'][1], i), lac(p['Ring_0'][2], p['Ring_Dihe'][2], i), lac(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh, oh, cz, ne1), lac(p['Ring_0'][0], p['Ring_Dihe'][0], i), lac(p['Ring_0'][1], p['Ring_Dihe'][1], i), lac(p['Ring_0'][2], p['Ring_Dihe'][2], i), lac(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(oh, cz, ce2, cz1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ne1, cz, ce2, cz1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz, ce2, cz1, hz1), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz, ce2, cz1, ch), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, ce3, cz2, hz2), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, ce3, cz2, ch), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, ce3, cz2, hz2), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(he3, ce3, cz2, ch), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce3, cz2, ch, hh2), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(ce3, cz2, ch, cz1), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz2, cz2, ch, hh2), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hz2, cz2, ch, cz1), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh2, ch, cz1, hz1), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hh2, ch, cz1, ce2), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz2, ch, cz1, hz1), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cz2, ch, cz1, ce2), lac(p['Ring_Dihe'][0], p['Ring_0'][0], i), lac(p['Ring_Dihe'][1], p['Ring_0'][1], i), lac(p['Ring_Dihe'][2], p['Ring_0'][2], i), lac(p['Ring_Dihe'][3], p['Ring_0'][3], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', ch, hh2), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', cz1, hz1), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', cz2, hz2), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', cz, oh), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', ce3, he3), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', ne1, he1), lac(p['Ami_imp'][0], p['Imp_0'][0], i), lac(p['Ami_imp'][1], p['Imp_0'][1], i), lac(p['Ami_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(cd2, 'CT', cg, cd1), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cg, lac(p['CA'][2], p['CA'][2], i), lac(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd1, lac(p['CA'][2], p['CA'][2], i), lac(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd1, lac(p['H4'][2], p['HA'][2], i), lac(p['H4'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ne1, lac(p['NA'][2], p['CA'][2], i), lac(p['NA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he1, lac(p['H'][2], p['HA'][2], i), lac(p['H'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ce2, lac(p['CA'][2], p['CA'][2], i), lac(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd2, lac(p['CA'][2], p['CA'][2], i), lac(p['CA'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ce3, lac(p['CA'][2], p['HA'][2], i), lac(p['CA'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), he3, lac(p['HA'][2], p['0_H'][2], i), lac(p['HA'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz, lac(p['0_C'][2], p['CA'][2], i), lac(p['0_C'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), oh, lac(p['0_O'][2], p['OH'][2], i), lac(p['0_O'][3], p['OH'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh, lac(p['0_H'][2], p['HO'][2], i), lac(p['0_H'][3], p['HO'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz1, lac(p['CA'][2], p['HA'][2], i), lac(p['CA'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz1, lac(p['HA'][2], p['0_H'][2], i), lac(p['HA'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cz2, lac(p['CA'][2], p['0_C'][2], i), lac(p['CA'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hz2, lac(p['HA'][2], p['0_H'][2], i), lac(p['HA'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ch, lac(p['CA'][2], p['0_C'][2], i), lac(p['CA'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hh2, lac(p['HA'][2], p['0_H'][2], i), lac(p['HA'][3], p['0_H'][3], i))
