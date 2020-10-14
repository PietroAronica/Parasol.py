# X3X to X4X Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/X3X.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/X4X.param', 'r') as b:
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
                change(parm, 'charge', ':{}@C3'.format(vxi), (fc['C3']/10)*i).execute()
                change(parm, 'charge', ':{}@H32'.format(vxi), (fc['H32']/10)*i).execute()
                change(parm, 'charge', ':{}@H33'.format(vxi), (fc['H33']/10)*i).execute()
                change(parm, 'charge', ':{}@C4'.format(vxi), bc['C3']+((fc['C4']-bc['C3'])/10)*i).execute()
                change(parm, 'charge', ':{}@H42'.format(vxi), bc['H32']+((fc['H42']-bc['H32'])/10)*i).execute()
                change(parm, 'charge', ':{}@H43'.format(vxi), bc['H33']+((fc['H43']-bc['H33'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        C2 = struct.residue_dict[aa].atom_dict['C2']
        C3 = struct.residue_dict[aa].atom_dict['C3']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'H23' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('C3', C2, C3))
                        	pdb.write(atom.superimposed1('H32', C3))
                        	pdb.write(atom.superimposed1('H33', C3))
			elif atom.get_name() == 'C3' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('C4'))
			elif atom.get_name() == 'H32' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('H42'))
			elif atom.get_name() == 'H33' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('H43'))
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
	ctrl.write("%s=loadpdb Param_files/LibPDB/X3X-X4X.pdb\n"%vxi)
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
	ctrl.write('set %s.1.1 type "CT"\n'%vxi)
	ctrl.write('set %s.1.2 type "HC"\n'%vxi)
	ctrl.write('set %s.1.3 type "HC"\n'%vxi)
	ctrl.write('set %s.1.4 type "CT"\n'%vxi)
	ctrl.write('set %s.1.5 type "HC"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, intcar))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, inthyd))
	ctrl.write('set %s.1.10 type "CT"\n'%vxi)
	ctrl.write('set %s.1.11 type "HC"\n'%vxi)
	ctrl.write('set %s.1.12 type "HC"\n'%vxi)
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
	ctrl.write('set %s.1 connect0 %s.1.C1\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C4\n'%(vxi, vxi))
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
