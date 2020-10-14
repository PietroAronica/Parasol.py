# XXX to XxX Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/XXX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/XxX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@CA1'.format(vxi), bc['CA1']+((fc['CA1']-bc['CA1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA12'.format(vxi), bc['HA12']+((fc['HA12']-bc['HA12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA13'.format(vxi), bc['HA13']+((fc['HA13']-bc['HA13'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB1'.format(vxi), bc['CB1']+((fc['CB1']-bc['CB1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB12'.format(vxi), bc['HB12']+((fc['HB12']-bc['HB12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB13'.format(vxi), bc['HB13']+((fc['HB13']-bc['HB13'])/10)*i).execute()
                change(parm, 'charge', ':{}@CC1'.format(vxi), (fc['CC1']/10)*i).execute()
                change(parm, 'charge', ':{}@HC12'.format(vxi), (fc['HC12']/10)*i).execute()
                change(parm, 'charge', ':{}@HC13'.format(vxi), (fc['HC13']/10)*i).execute()
                change(parm, 'charge', ':{}@C1'.format(vxi), bc['C1']+((fc['C1']-bc['C1'])/10)*i).execute()
                change(parm, 'charge', ':{}@C21'.format(vxi), bc['C21']+((fc['C21']-bc['C21'])/10)*i).execute()
                change(parm, 'charge', ':{}@H21'.format(vxi), bc['H21']+((fc['H21']-bc['H21'])/10)*i).execute()
                change(parm, 'charge', ':{}@C31'.format(vxi), bc['C31']+((fc['C31']-bc['C31'])/10)*i).execute()
                change(parm, 'charge', ':{}@H31'.format(vxi), bc['H31']+((fc['H31']-bc['H31'])/10)*i).execute()
                change(parm, 'charge', ':{}@C4'.format(vxi), bc['C4']+((fc['C4']-bc['C4'])/10)*i).execute()
                change(parm, 'charge', ':{}@C32'.format(vxi), bc['C32']+((fc['C32']-bc['C32'])/10)*i).execute()
                change(parm, 'charge', ':{}@H32'.format(vxi), bc['H32']+((fc['H32']-bc['H32'])/10)*i).execute()
                change(parm, 'charge', ':{}@C22'.format(vxi), bc['C22']+((fc['C22']-bc['C22'])/10)*i).execute()
                change(parm, 'charge', ':{}@H22'.format(vxi), bc['H22']+((fc['H22']-bc['H22'])/10)*i).execute()
                change(parm, 'charge', ':{}@CC2'.format(vxi), (fc['CC2']/10)*i).execute()
                change(parm, 'charge', ':{}@HC22'.format(vxi), (fc['HC22']/10)*i).execute()
                change(parm, 'charge', ':{}@HC23'.format(vxi), (fc['HC23']/10)*i).execute()
                change(parm, 'charge', ':{}@CB2'.format(vxi), bc['CB2']+((fc['CB2']-bc['CB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB22'.format(vxi), bc['HB22']+((fc['HB22']-bc['HB22'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB23'.format(vxi), bc['HB23']+((fc['HB23']-bc['HB23'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA2'.format(vxi), bc['CA2']+((fc['CA2']-bc['CA2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA22'.format(vxi), bc['HA22']+((fc['HA22']-bc['HA22'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA23'.format(vxi), bc['HA23']+((fc['HA23']-bc['HA23'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CB1 = struct.residue_dict[aa].atom_dict['CB1']
        C1 = struct.residue_dict[aa].atom_dict['C1']
        CB2 = struct.residue_dict[aa].atom_dict['CB2']
        C4 = struct.residue_dict[aa].atom_dict['C4']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HB13' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('CC1', CB1, C1))
                        	pdb.write(atom.superimposed1('HC12', CB1))
                        	pdb.write(atom.superimposed1('HC13', CB1))
			elif atom.get_name() == 'H22' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('CC2', CB2, C4))
                        	pdb.write(atom.superimposed1('HC22', CB2))
                        	pdb.write(atom.superimposed1('HC23', CB2))
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
	return var1, var2, var3

def lib_make(ff, outputfile, vxi='VXI', var=variablemake()):
 	intcar = var[0]
 	inthyd = var[1]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/XXX-XxX.pdb\n"%vxi)
	ctrl.write('set %s.1.1  element "C"\n'%vxi)
	ctrl.write('set %s.1.2  element "H"\n'%vxi)
	ctrl.write('set %s.1.3  element "H"\n'%vxi)
	ctrl.write('set %s.1.4  element "C"\n'%vxi)
	ctrl.write('set %s.1.5  element "H"\n'%vxi)
	ctrl.write('set %s.1.6  element "H"\n'%vxi)
	ctrl.write('set %s.1.7  element "C"\n'%vxi)
	ctrl.write('set %s.1.8  element "H"\n'%vxi)
	ctrl.write('set %s.1.9  element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "C"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "C"\n'%vxi)
	ctrl.write('set %s.1.24 element "H"\n'%vxi)
	ctrl.write('set %s.1.25 element "H"\n'%vxi)
	ctrl.write('set %s.1.26 element "C"\n'%vxi)
	ctrl.write('set %s.1.27 element "H"\n'%vxi)
	ctrl.write('set %s.1.28 element "H"\n'%vxi)
	ctrl.write('set %s.1.1  name "CA1"\n'%vxi)
	ctrl.write('set %s.1.2  name "HA12"\n'%vxi)
	ctrl.write('set %s.1.3  name "HA13"\n'%vxi)
	ctrl.write('set %s.1.4  name "CB1"\n'%vxi)
	ctrl.write('set %s.1.5  name "HB12"\n'%vxi)
	ctrl.write('set %s.1.6  name "HB13"\n'%vxi)
	ctrl.write('set %s.1.7  name "CC1"\n'%vxi)
	ctrl.write('set %s.1.8  name "HC12"\n'%vxi)
	ctrl.write('set %s.1.9  name "HC13"\n'%vxi)
	ctrl.write('set %s.1.10 name "C1"\n'%vxi)
	ctrl.write('set %s.1.11 name "C21"\n'%vxi)
	ctrl.write('set %s.1.12 name "H21"\n'%vxi)
	ctrl.write('set %s.1.13 name "C31"\n'%vxi)
	ctrl.write('set %s.1.14 name "H31"\n'%vxi)
	ctrl.write('set %s.1.15 name "C4"\n'%vxi)
	ctrl.write('set %s.1.16 name "C32"\n'%vxi)
	ctrl.write('set %s.1.17 name "H32"\n'%vxi)
	ctrl.write('set %s.1.18 name "C22"\n'%vxi)
	ctrl.write('set %s.1.19 name "H22"\n'%vxi)
	ctrl.write('set %s.1.20 name "CC2"\n'%vxi)
	ctrl.write('set %s.1.21 name "HC22"\n'%vxi)
	ctrl.write('set %s.1.22 name "HC23"\n'%vxi)
	ctrl.write('set %s.1.23 name "CB2"\n'%vxi)
	ctrl.write('set %s.1.24 name "HB22"\n'%vxi)
	ctrl.write('set %s.1.25 name "HB23"\n'%vxi)
	ctrl.write('set %s.1.26 name "CA2"\n'%vxi)
	ctrl.write('set %s.1.27 name "HA22"\n'%vxi)
	ctrl.write('set %s.1.28 name "HA23"\n'%vxi)
	ctrl.write('set %s.1.1  type "CT"\n'%vxi)
	ctrl.write('set %s.1.2  type "HC"\n'%vxi)
	ctrl.write('set %s.1.3  type "HC"\n'%vxi)
	ctrl.write('set %s.1.4  type "CT"\n'%vxi)
	ctrl.write('set %s.1.5  type "HC"\n'%vxi)
	ctrl.write('set %s.1.6  type "HC"\n'%vxi)
	ctrl.write('set %s.1.7  type "%s"\n'%(vxi, intcar))
	ctrl.write('set %s.1.8  type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.9  type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.10 type "CA"\n'%vxi)
	ctrl.write('set %s.1.11 type "CA"\n'%vxi)
	ctrl.write('set %s.1.12 type "HA"\n'%vxi)
	ctrl.write('set %s.1.13 type "CA"\n'%vxi)
	ctrl.write('set %s.1.14 type "HA"\n'%vxi)
	ctrl.write('set %s.1.15 type "CA"\n'%vxi)
	ctrl.write('set %s.1.16 type "CA"\n'%vxi)
	ctrl.write('set %s.1.17 type "HA"\n'%vxi)
	ctrl.write('set %s.1.18 type "CA"\n'%vxi)
	ctrl.write('set %s.1.19 type "HA"\n'%vxi)
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, intcar))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.23 type "CT"\n'%vxi)
	ctrl.write('set %s.1.24 type "HC"\n'%vxi)
	ctrl.write('set %s.1.25 type "HC"\n'%vxi)
	ctrl.write('set %s.1.26 type "CT"\n'%vxi)
	ctrl.write('set %s.1.27 type "HC"\n'%vxi)
	ctrl.write('set %s.1.28 type "HC"\n'%vxi)
	ctrl.write('bond %s.1.1  %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1  %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1  %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4  %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4  %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4  %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7  %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7  %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7  %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.26\n'%(vxi, vxi))
	ctrl.write('bond %s.1.26 %s.1.27\n'%(vxi, vxi))
	ctrl.write('bond %s.1.26 %s.1.28\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.CA1\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.CA2\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.CA1\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.CA2\n'%(vxi, vxi))
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
	intcar = var[0]
	inthyd = var[1]
        Frcmod_creator.make_hyb()
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), intcar, cal(p['0_C'][0], p['CT'][0], i), cal(p['0_C'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), inthyd, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(intcar, inthyd), cal(p['CT_mH'][0], p['CT_HC'][0], i), cal(p['CT_mH'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(intcar, 'CT'), cal(p['CT_mCT'][0], p['CT_CT'][0], i), cal(p['CT_mCT'][1], p['CT_CT'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(intcar, 'CA'), cal(p['CT_mCT'][0], p['CT_CA'][0], i), cal(p['CT_mCT'][1], p['CT_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', 'CT', intcar), cal(p['C_C_S'][0], p['C_C_S'][0], i), cal(p['C_C_S'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', intcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', intcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', intcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intcar, inthyd), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intcar, 'CA'), cal(p['Dritt'][0], p['C_C_CA'][0], i), cal(p['Dritt'][1], p['C_C_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(inthyd, intcar, 'CA'), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(inthyd, intcar, inthyd), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CA', 'CA', intcar), cal(p['F_CT_CA_CA'][0], p['F_CT_CA_CA'][0], i), cal(p['F_CT_CA_CA'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(inthyd, intcar, 'CT', 'HC'), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(inthyd, intcar, 'CT', 'H1'), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(inthyd, intcar, 'CT', 'CT'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intcar, 'CT', 'H1'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', intcar, inthyd), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', intcar, 'CT'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', intcar, 'CA'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', intcar, 'CA'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', intcar, 'CA'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', intcar, 'CA'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', intcar, 'CT'), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', intcar, inthyd), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CA', 'CA', 'CA', intcar), cal(p['Imp_0'][0], p['Ring_imp'][0], i), cal(p['Imp_0'][1], p['Ring_imp'][1], i), cal(p['Imp_0'][2], p['Ring_imp'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), intcar, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), inthyd, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
