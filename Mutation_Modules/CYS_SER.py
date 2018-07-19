# CYS to SER Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/CYS.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/SER.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@OG'.format(vxi), bc['SG']+((fc['OG']-bc['SG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG'.format(vxi), bc['HG']+((fc['HG']-bc['HG'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'SG' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('OG'))
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
	thioxy = var[0]
	thihyd = var[1]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/SER-CYS.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "O"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "OG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG"\n'%vxi)
	ctrl.write('set %s.1.10 name "C"\n'%vxi)
	ctrl.write('set %s.1.11 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "H1"\n'%vxi)
	ctrl.write('set %s.1.7 type "H1"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, thioxy))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, thihyd))
	ctrl.write('set %s.1.10 type "C"\n'%vxi)
	ctrl.write('set %s.1.11 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
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

def stock_add_to_all(var=variablemake()):
	thioxy = var[0]
	thihyd = var[1]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(thioxy, 'O', 'sp3')
	Frcmod_creator.TYPE_insert(thihyd, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), thioxy, lac(p['OH'][0], p['SH'][0], i), lac(p['OH'][1], p['SH'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), thihyd, lac(p['HO'][0], p['HS'][0], i), lac(p['HO'][1], p['HS'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', thioxy), lac(p['CT_OH'][0], p['CT_SH'][0], i), lac(p['CT_OH'][1], p['CT_SH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(thioxy, thihyd), lac(p['OH_HO'][0], p['SH_HS'][0], i), lac(p['OH_HO'][1], p['SH_HS'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', thioxy), lac(p['C_C_H'][0], p['C_C_SH'][0], i), lac(p['C_C_H'][1], p['C_C_SH'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', thioxy), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', thioxy, thihyd), lac(p['C_O_H'][0], p['C_SH_H'][0], i), lac(p['C_O_H'][1], p['C_SH_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(thihyd, thioxy, 'CT', 'H1'), lac(p['X_C_O_X'][0], p['X_C_SH_X'][0], i), lac(p['X_C_O_X'][1], p['X_C_SH_X'][1], i), lac(p['X_C_O_X'][2], p['X_C_SH_X'][2], i), lac(p['X_C_O_X'][3], p['X_C_SH_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(thihyd, thioxy, 'CT', 'CT'), lac(p['C_C_O_H_2'][0], p['0_3'][0], i), lac(p['C_C_O_H_2'][1], p['0_3'][1], i), lac(p['C_C_O_H_2'][2], p['0_3'][2], i), lac(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(thihyd, thioxy, 'CT', 'CT'), lac(p['C_C_O_H_1'][0], p['X_C_SH_X'][0], i), lac(p['C_C_O_H_1'][1], p['X_C_SH_X'][1], i), lac(p['C_C_O_H_1'][2], p['X_C_SH_X'][2], i), lac(p['C_C_O_H_1'][3], p['X_C_SH_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(thioxy, 'CT', 'CT', 'H1'), lac(p['C_C_O_H_2'][0], p['0_3'][0], i), lac(p['C_C_O_H_2'][1], p['0_3'][1], i), lac(p['C_C_O_H_2'][2], p['0_3'][2], i), lac(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(thioxy, 'CT', 'CT', 'H1'), lac(p['C_C_O_H_1'][0], p['X_C_C_X'][0], i), lac(p['C_C_O_H_1'][1], p['X_C_C_X'][1], i), lac(p['C_C_O_H_1'][2], p['X_C_C_X'][2], i), lac(p['C_C_O_H_1'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), thioxy, lac(p['OH'][2], p['SH'][2], i), lac(p['OH'][3], p['SH'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), thihyd, lac(p['HO'][2], p['HS'][2], i), lac(p['HO'][3], p['HS'][3], i))
