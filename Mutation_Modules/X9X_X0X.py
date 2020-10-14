# X9X to X0X Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/X9X.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/X0X.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@C1'.format(vxi), bc['C1']+((fc['C1']-bc['C1'])/10)*i).execute()
                change(parm, 'charge', ':{}@H12'.format(vxi), bc['H12']+((fc['H12']-bc['H12'])/10)*i).execute()
                change(parm, 'charge', ':{}@H13'.format(vxi), bc['H13']+((fc['H13']-bc['H13'])/10)*i).execute()
                change(parm, 'charge', ':{}@C2'.format(vxi), bc['C2']+((fc['C2']-bc['C2'])/10)*i).execute()
                change(parm, 'charge', ':{}@H22'.format(vxi), bc['H22']+((fc['H22']-bc['H22'])/10)*i).execute()
                change(parm, 'charge', ':{}@H23'.format(vxi), bc['H23']+((fc['H23']-bc['H23'])/10)*i).execute()
                change(parm, 'charge', ':{}@C3'.format(vxi), bc['C3']+((fc['C3']-bc['C3'])/10)*i).execute()
                change(parm, 'charge', ':{}@H32'.format(vxi), bc['H32']+((fc['H32']-bc['H32'])/10)*i).execute()
                change(parm, 'charge', ':{}@H33'.format(vxi), bc['H33']+((fc['H33']-bc['H33'])/10)*i).execute()
                change(parm, 'charge', ':{}@C4'.format(vxi), bc['C4']+((fc['C4']-bc['C4'])/10)*i).execute()
                change(parm, 'charge', ':{}@H42'.format(vxi), bc['H42']+((fc['H42']-bc['H42'])/10)*i).execute()
                change(parm, 'charge', ':{}@H43'.format(vxi), bc['H43']+((fc['H43']-bc['H43'])/10)*i).execute()
                change(parm, 'charge', ':{}@C5'.format(vxi), bc['C5']+((fc['C5']-bc['C5'])/10)*i).execute()
                change(parm, 'charge', ':{}@H52'.format(vxi), bc['H52']+((fc['H52']-bc['H52'])/10)*i).execute()
                change(parm, 'charge', ':{}@H53'.format(vxi), bc['H53']+((fc['H53']-bc['H53'])/10)*i).execute()
                change(parm, 'charge', ':{}@C6'.format(vxi), bc['C6']+((fc['C6']-bc['C6'])/10)*i).execute()
                change(parm, 'charge', ':{}@H62'.format(vxi), bc['H62']+((fc['H62']-bc['H62'])/10)*i).execute()
                change(parm, 'charge', ':{}@H63'.format(vxi), bc['H63']+((fc['H63']-bc['H63'])/10)*i).execute()
                change(parm, 'charge', ':{}@C7'.format(vxi), bc['C7']+((fc['C7']-bc['C7'])/10)*i).execute()
                change(parm, 'charge', ':{}@H72'.format(vxi), bc['H72']+((fc['H72']-bc['H72'])/10)*i).execute()
                change(parm, 'charge', ':{}@H73'.format(vxi), bc['H73']+((fc['H73']-bc['H73'])/10)*i).execute()
                change(parm, 'charge', ':{}@C8'.format(vxi), bc['C8']+((fc['C8']-bc['C8'])/10)*i).execute()
                change(parm, 'charge', ':{}@H82'.format(vxi), bc['H82']+((fc['H82']-bc['H82'])/10)*i).execute()
                change(parm, 'charge', ':{}@H83'.format(vxi), bc['H83']+((fc['H83']-bc['H83'])/10)*i).execute()
                change(parm, 'charge', ':{}@C9'.format(vxi), (fc['C9']/10)*i).execute()
                change(parm, 'charge', ':{}@H92'.format(vxi), (fc['H92']/10)*i).execute()
                change(parm, 'charge', ':{}@H93'.format(vxi), (fc['H93']/10)*i).execute()
                change(parm, 'charge', ':{}@C10'.format(vxi), bc['C9']+((fc['C10']-bc['C9'])/10)*i).execute()
                change(parm, 'charge', ':{}@H102'.format(vxi), bc['H92']+((fc['H102']-bc['H92'])/10)*i).execute()
                change(parm, 'charge', ':{}@H103'.format(vxi), bc['H93']+((fc['H103']-bc['H93'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        C8 = struct.residue_dict[aa].atom_dict['C8']
        C9 = struct.residue_dict[aa].atom_dict['C9']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'H83' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('C9', C8, C9))
                        	pdb.write(atom.superimposed1('H92', C9))
                        	pdb.write(atom.superimposed2('H93', C9))
			elif atom.get_name() == 'C9' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('C10'))
			elif atom.get_name() == 'H92' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('H102'))
			elif atom.get_name() == 'H93' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('H103'))
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
	ctrl.write("%s=loadpdb Param_files/LibPDB/X9X-X0X.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "C"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "H"\n'%vxi)
	ctrl.write('set %s.1.4 element "C"\n'%vxi)
	ctrl.write('set %s.1.5 element "H"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "C"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "C"\n'%vxi)
	ctrl.write('set %s.1.23 element "H"\n'%vxi)
	ctrl.write('set %s.1.24 element "H"\n'%vxi)
	ctrl.write('set %s.1.25 element "C"\n'%vxi)
	ctrl.write('set %s.1.26 element "H"\n'%vxi)
	ctrl.write('set %s.1.27 element "H"\n'%vxi)
	ctrl.write('set %s.1.28 element "C"\n'%vxi)
	ctrl.write('set %s.1.29 element "H"\n'%vxi)
	ctrl.write('set %s.1.30 element "H"\n'%vxi)
	ctrl.write('set %s.1.1 name "C1"\n'%vxi)
	ctrl.write('set %s.1.2 name "H12"\n'%vxi)
	ctrl.write('set %s.1.3 name "H13"\n'%vxi)
	ctrl.write('set %s.1.4 name "C2"\n'%vxi)
	ctrl.write('set %s.1.5 name "H22"\n'%vxi)
	ctrl.write('set %s.1.6 name "H23"\n'%vxi)
	ctrl.write('set %s.1.7 name "C3"\n'%vxi)
	ctrl.write('set %s.1.8 name "H32"\n'%vxi)
	ctrl.write('set %s.1.9 name "H33"\n'%vxi)
	ctrl.write('set %s.1.10 name "C4"\n'%vxi)
	ctrl.write('set %s.1.11 name "H42"\n'%vxi)
	ctrl.write('set %s.1.12 name "H43"\n'%vxi)
	ctrl.write('set %s.1.13 name "C5"\n'%vxi)
	ctrl.write('set %s.1.14 name "H52"\n'%vxi)
	ctrl.write('set %s.1.15 name "H53"\n'%vxi)
	ctrl.write('set %s.1.16 name "C6"\n'%vxi)
	ctrl.write('set %s.1.17 name "H62"\n'%vxi)
	ctrl.write('set %s.1.18 name "H63"\n'%vxi)
	ctrl.write('set %s.1.19 name "C7"\n'%vxi)
	ctrl.write('set %s.1.20 name "H72"\n'%vxi)
	ctrl.write('set %s.1.21 name "H73"\n'%vxi)
	ctrl.write('set %s.1.22 name "C8"\n'%vxi)
	ctrl.write('set %s.1.23 name "H82"\n'%vxi)
	ctrl.write('set %s.1.24 name "H83"\n'%vxi)
	ctrl.write('set %s.1.25 name "C9"\n'%vxi)
	ctrl.write('set %s.1.26 name "H92"\n'%vxi)
	ctrl.write('set %s.1.27 name "H93"\n'%vxi)
	ctrl.write('set %s.1.28 name "C10"\n'%vxi)
	ctrl.write('set %s.1.29 name "H102"\n'%vxi)
	ctrl.write('set %s.1.30 name "H103"\n'%vxi)
	ctrl.write('set %s.1.1 type "CT"\n'%vxi)
	ctrl.write('set %s.1.2 type "HC"\n'%vxi)
	ctrl.write('set %s.1.3 type "HC"\n'%vxi)
	ctrl.write('set %s.1.4 type "CT"\n'%vxi)
	ctrl.write('set %s.1.5 type "HC"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "CT"\n'%vxi)
	ctrl.write('set %s.1.8 type "HC"\n'%vxi)
	ctrl.write('set %s.1.9 type "HC"\n'%vxi)
	ctrl.write('set %s.1.10 type "CT"\n'%vxi)
	ctrl.write('set %s.1.11 type "HC"\n'%vxi)
	ctrl.write('set %s.1.12 type "HC"\n'%vxi)
	ctrl.write('set %s.1.13 type "CT"\n'%vxi)
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "HC"\n'%vxi)
	ctrl.write('set %s.1.16 type "CT"\n'%vxi)
	ctrl.write('set %s.1.17 type "HC"\n'%vxi)
	ctrl.write('set %s.1.18 type "HC"\n'%vxi)
	ctrl.write('set %s.1.19 type "CT"\n'%vxi)
	ctrl.write('set %s.1.20 type "HC"\n'%vxi)
	ctrl.write('set %s.1.21 type "HC"\n'%vxi)
	ctrl.write('set %s.1.22 type "CT"\n'%vxi)
	ctrl.write('set %s.1.23 type "HC"\n'%vxi)
	ctrl.write('set %s.1.24 type "HC"\n'%vxi)
	ctrl.write('set %s.1.25 type "%s"\n'%(vxi, intcar))
	ctrl.write('set %s.1.26 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.27 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.28 type "CT"\n'%vxi)
	ctrl.write('set %s.1.29 type "HC"\n'%vxi)
	ctrl.write('set %s.1.30 type "HC"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.22 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.22 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.22 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.26\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.27\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.28\n'%(vxi, vxi))
	ctrl.write('bond %s.1.28 %s.1.29\n'%(vxi, vxi))
	ctrl.write('bond %s.1.28 %s.1.30\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.C1\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C10\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.C1\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C10\n'%(vxi, vxi))
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
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', 'CT', intcar), cal(p['C_C_S'][0], p['C_C_S'][0], i), cal(p['C_C_S'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', intcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', intcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', intcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intcar, inthyd), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', intcar, 'CT'), cal(p['Dritt'][0], p['C_C_C'][0], i), cal(p['Dritt'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(inthyd, intcar, inthyd), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Close'][1], p['H_C_H'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(inthyd, intcar, 'CT', 'HC'), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(inthyd, intcar, 'CT', 'H1'), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(inthyd, intcar, 'CT', 'CT'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intcar, 'CT', 'H1'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intcar, 'CT', 'HC'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', intcar, inthyd), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', intcar, 'CT'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intcar, 'CT', 'CT'), cal(p['0_12'][0], p['C_C_C_C_1'][0], i), cal(p['0_12'][1], p['C_C_C_C_1'][1], i), cal(p['0_12'][2], p['C_C_C_C_1'][2], i), cal(p['0_12'][3], p['C_C_C_C_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intcar, 'CT', 'CT'), cal(p['0_11'][0], p['C_C_C_C_2'][0], i), cal(p['0_11'][1], p['C_C_C_C_2'][1], i), cal(p['0_11'][2], p['C_C_C_C_2'][2], i), cal(p['0_11'][3], p['C_C_C_C_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', intcar, 'CT', 'CT'), cal(p['0_2'][0], p['C_C_C_C_3'][0], i), cal(p['0_2'][1], p['C_C_C_C_3'][1], i), cal(p['0_2'][2], p['C_C_C_C_3'][2], i), cal(p['0_2'][3], p['C_C_C_C_3'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), intcar, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), inthyd, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
