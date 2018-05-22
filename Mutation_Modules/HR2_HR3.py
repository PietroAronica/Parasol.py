# HR2 to HR3 Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/HR2.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/HR3.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		i = float(i)
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
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
                change(parm, 'charge', ':{}@CE'.format(vxi), bc['CE']+((fc['CE']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), bc['HE2']+((fc['HE2']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), bc['HE3']+((fc['HE3']-bc['HE3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ'.format(vxi), fc['CZ']/10*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vxi), fc['HZ2']/10*i).execute()
                change(parm, 'charge', ':{}@HZ3'.format(vxi), fc['HZ3']/10*i).execute()
                change(parm, 'charge', ':{}@CH'.format(vxi), bc['CZ']+((fc['CH']-bc['CZ'])/10)*i).execute()
                change(parm, 'charge', ':{}@HH2'.format(vxi), bc['HZ2']+((fc['HH2']-bc['HZ2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HH3'.format(vxi), bc['HZ3']+((fc['HH3']-bc['HZ3'])/10)*i).execute()
                change(parm, 'charge', ':{}@NT'.format(vxi), bc['NH']+((fc['NT']-bc['NH'])/10)*i).execute()
                change(parm, 'charge', ':{}@HT'.format(vxi), bc['HH']+((fc['HT']-bc['HH'])/10)*i).execute()
                change(parm, 'charge', ':{}@CI'.format(vxi), bc['CT']+((fc['CI']-bc['CT'])/10)*i).execute()
                change(parm, 'charge', ':{}@NK1'.format(vxi), bc['NI1']+((fc['NK1']-bc['NI1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HK11'.format(vxi), bc['HI11']+((fc['HK11']-bc['HI11'])/10)*i).execute()
                change(parm, 'charge', ':{}@HK12'.format(vxi), bc['HI12']+((fc['HK12']-bc['HI12'])/10)*i).execute()
                change(parm, 'charge', ':{}@NK2'.format(vxi), bc['NI2']+((fc['NK2']-bc['NI2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HK21'.format(vxi), bc['HI21']+((fc['HK21']-bc['HI21'])/10)*i).execute()
                change(parm, 'charge', ':{}@HK22'.format(vxi), bc['HI22']+((fc['HK22']-bc['HI22'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		#print printDetails(parm, ':VXI')
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CE = struct.residue_dict[aa].atom_dict['CE']
        CZ = struct.residue_dict[aa].atom_dict['CZ']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HE3' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('CZ', CE, CZ))
                        	pdb.write(atom.superimposed1('HZ2', CZ))
                        	pdb.write(atom.superimposed2('HZ3', CZ))
			elif atom.get_name() == 'CZ' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CH'))
			elif atom.get_name() == 'HZ2' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HH2'))
			elif atom.get_name() == 'HZ3' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HH3'))
			elif atom.get_name() == 'NH' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('NT'))
			elif atom.get_name() == 'HH' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HT'))
			elif atom.get_name() == 'CT' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CI'))
			elif atom.get_name() == 'NI1' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('NK1'))
			elif atom.get_name() == 'HI11' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HK11'))
			elif atom.get_name() == 'HI12' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HK12'))
			elif atom.get_name() == 'NI2' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('NK2'))
			elif atom.get_name() == 'HI21' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HK21'))
			elif atom.get_name() == 'HI22' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HK22'))
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

def lib_make(ff, outputfile, vxi='VXI', intcar='ic', inthyd='ih', norcar='ct'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/HR2-HR3.pdb\n"%vxi)
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
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "C"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "C"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "N"\n'%vxi)
	ctrl.write('set %s.1.24 element "H"\n'%vxi)
	ctrl.write('set %s.1.25 element "C"\n'%vxi)
	ctrl.write('set %s.1.26 element "N"\n'%vxi)
	ctrl.write('set %s.1.27 element "H"\n'%vxi)
	ctrl.write('set %s.1.28 element "H"\n'%vxi)
	ctrl.write('set %s.1.29 element "N"\n'%vxi)
	ctrl.write('set %s.1.30 element "H"\n'%vxi)
	ctrl.write('set %s.1.31 element "H"\n'%vxi)
	ctrl.write('set %s.1.32 element "C"\n'%vxi)
	ctrl.write('set %s.1.33 element "O"\n'%vxi)
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
	ctrl.write('set %s.1.14 name "CE"\n'%vxi)
	ctrl.write('set %s.1.15 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.16 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.17 name "CZ"\n'%vxi)
	ctrl.write('set %s.1.18 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.19 name "HZ3"\n'%vxi)
	ctrl.write('set %s.1.20 name "CH"\n'%vxi)
	ctrl.write('set %s.1.21 name "HH2"\n'%vxi)
	ctrl.write('set %s.1.22 name "HH3"\n'%vxi)
	ctrl.write('set %s.1.23 name "NT"\n'%vxi)
	ctrl.write('set %s.1.24 name "HT"\n'%vxi)
	ctrl.write('set %s.1.25 name "CI"\n'%vxi)
	ctrl.write('set %s.1.26 name "NK1"\n'%vxi)
	ctrl.write('set %s.1.27 name "HK11"\n'%vxi)
	ctrl.write('set %s.1.28 name "HK12"\n'%vxi)
	ctrl.write('set %s.1.29 name "NK2"\n'%vxi)
	ctrl.write('set %s.1.30 name "HK21"\n'%vxi)
	ctrl.write('set %s.1.31 name "HK22"\n'%vxi)
	ctrl.write('set %s.1.32 name "C"\n'%vxi)
	ctrl.write('set %s.1.33 name "O"\n'%vxi)
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
	ctrl.write('set %s.1.12 type "HC"\n'%vxi)
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, norcar))
	ctrl.write('set %s.1.15 type "HC"\n'%vxi)
	ctrl.write('set %s.1.16 type "HC"\n'%vxi)
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, intcar))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.20 type "CT"\n'%vxi)
	ctrl.write('set %s.1.21 type "H1"\n'%vxi)
	ctrl.write('set %s.1.22 type "H1"\n'%vxi)
	ctrl.write('set %s.1.23 type "N2"\n'%vxi)
	ctrl.write('set %s.1.24 type "H"\n'%vxi)
	ctrl.write('set %s.1.25 type "CA"\n'%vxi)
	ctrl.write('set %s.1.26 type "N2"\n'%vxi)
	ctrl.write('set %s.1.27 type "H"\n'%vxi)
	ctrl.write('set %s.1.28 type "H"\n'%vxi)
	ctrl.write('set %s.1.29 type "N2"\n'%vxi)
	ctrl.write('set %s.1.30 type "H"\n'%vxi)
	ctrl.write('set %s.1.31 type "H"\n'%vxi)
	ctrl.write('set %s.1.32 type "C"\n'%vxi)
	ctrl.write('set %s.1.33 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.32\n'%(vxi, vxi))
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
	ctrl.write('bond %s.1.20 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.26\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.29\n'%(vxi, vxi))
	ctrl.write('bond %s.1.26 %s.1.27\n'%(vxi, vxi))
	ctrl.write('bond %s.1.26 %s.1.28\n'%(vxi, vxi))
	ctrl.write('bond %s.1.29 %s.1.30\n'%(vxi, vxi))
	ctrl.write('bond %s.1.29 %s.1.31\n'%(vxi, vxi))
	ctrl.write('bond %s.1.32 %s.1.33\n'%(vxi, vxi))
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

def lac(y, x, i):
	num = x+((y-x)/10)*i
	return num

def stock_add_to_all(intcar='ic', inthyd='ih', norcar='ct'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(norcar, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(intcar, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(inthyd, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), norcar, cal(p['CT'][0], p['CT'][0], i), cal(p['CT'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), intcar, cal(p['0_C'][0], p['CT'][0], i), cal(p['0_C'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), inthyd, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(norcar, intcar), cal(p['CT_mC'][0], p['CT_CT'][0], i), cal(p['CT_mC'][1], p['CT_CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', intcar), cal(p['CT_mC'][0], p['CT_CT'][0], i), cal(p['CT_mC'][1], p['CT_CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', norcar), cal(p['CT_CT'][0], p['CT_CT'][0], i), cal(p['CT_CT'][1], p['CT_CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(norcar, 'HC'), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(intcar, inthyd), cal(p['HC_mC'][0], p['CT_HC'][0], i), cal(p['HC_mC'][1], p['CT_HC'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(inthyd, intcar, inthyd), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(norcar, intcar, inthyd), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intcar, inthyd), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(norcar, intcar, 'CT'), cal(p['Dritt'][0], p['C_C_C'][0], i), cal(p['Dritt'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', norcar, intcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', norcar, intcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', norcar, 'CT'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', norcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(norcar, 'CT', 'CT'), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', norcar, 'HC'), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intcar, 'CT', 'H1'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(intcar, 'CT', 'N2'), cal(p['C_C_N2'][0], p['C_C_N2'][0], i), cal(p['C_C_N2'][1], p['C_C_N2'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', intcar, inthyd), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N2', 'CT', intcar, inthyd), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', norcar, intcar, inthyd), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', norcar, intcar, inthyd), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', intcar, norcar), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N2', 'CT', intcar, norcar), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', norcar, intcar, 'CT'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', norcar, intcar, 'CT'), cal(p['0_12'][0], p['C_C_C_C_1'][0], i), cal(p['0_12'][1], p['C_C_C_C_1'][1], i), cal(p['0_12'][2], p['C_C_C_C_1'][2], i), cal(p['0_12'][3], p['C_C_C_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', norcar, intcar, 'CT'), cal(p['0_11'][0], p['C_C_C_C_2'][0], i), cal(p['0_11'][1], p['C_C_C_C_2'][1], i), cal(p['0_11'][2], p['C_C_C_C_2'][2], i), cal(p['0_11'][3], p['C_C_C_C_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', norcar, intcar, 'CT'), cal(p['0_2'][0], p['C_C_C_C_3'][0], i), cal(p['0_2'][1], p['C_C_C_C_3'][1], i), cal(p['0_2'][2], p['C_C_C_C_3'][2], i), cal(p['0_2'][3], p['C_C_C_C_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', norcar, intcar), cal(p['0_12'][0], p['C_C_C_C_1'][0], i), cal(p['0_12'][1], p['C_C_C_C_1'][1], i), cal(p['0_12'][2], p['C_C_C_C_1'][2], i), cal(p['0_12'][3], p['C_C_C_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', norcar, intcar), cal(p['0_11'][0], p['C_C_C_C_2'][0], i), cal(p['0_11'][1], p['C_C_C_C_2'][1], i), cal(p['0_11'][2], p['C_C_C_C_2'][2], i), cal(p['0_11'][3], p['C_C_C_C_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', norcar, intcar), cal(p['0_2'][0], p['C_C_C_C_3'][0], i), cal(p['0_2'][1], p['C_C_C_C_3'][1], i), cal(p['0_2'][2], p['C_C_C_C_3'][2], i), cal(p['0_2'][3], p['C_C_C_C_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', norcar, 'HC'), cal(p['H_C_C_H'][0], p['H_C_C_H'][0], i), cal(p['H_C_C_H'][1], p['H_C_C_H'][1], i), cal(p['H_C_C_H'][2], p['H_C_C_H'][2], i), cal(p['H_C_C_H'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', norcar, 'HC'), cal(p['C_C_C_H'][0], p['C_C_C_H'][0], i), cal(p['C_C_C_H'][1], p['C_C_C_H'][1], i), cal(p['C_C_C_H'][2], p['C_C_C_H'][2], i), cal(p['C_C_C_H'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', norcar, intcar), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), norcar, cal(p['CT'][2], p['CT'][2], i), cal(p['CT'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), intcar, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), inthyd, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
