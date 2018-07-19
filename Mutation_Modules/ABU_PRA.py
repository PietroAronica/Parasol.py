# ABU to PRA Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/ABU.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/PRA.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HG1 :{}@HD 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HG2 :{}@HD 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vxi), bc['HG1']-(bc['HG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG2']-(bc['HG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HG2']-(bc['HG3']/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), (fc['CD']/10)*i).execute()
                change(parm, 'charge', ':{}@HD'.format(vxi), (fc['HD']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CG = struct.residue_dict[aa].atom_dict['CG']
        HG1 = struct.residue_dict[aa].atom_dict['HG1']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HG3' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('CD', CG, HG1))
                        	pdb.write(atom.superimposed1('HD', HG1))
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
	return var1, var2, var3, var4, var5, var6, var7, var8

def lib_make(ff, outputfile, vxi='VXI', var=variablemake()):
	difcar = var[0]
	hyd1 = var[1]
	hyd2 = var[2]
	newcar = var[3]
	newhyd = var[4]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ABU-PRA.pdb\n"%vxi)
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
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.12 name "CD"\n'%vxi)
	ctrl.write('set %s.1.13 name "HD"\n'%vxi)
	ctrl.write('set %s.1.14 name "C"\n'%vxi)
	ctrl.write('set %s.1.15 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, difcar))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, hyd1))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, hyd2))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, hyd2))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.14 type "C"\n'%vxi)
	ctrl.write('set %s.1.15 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
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
	difcar = var[0]
	hyd1 = var[1]
	hyd2 = var[2]
	newcar = var[3]
	newhyd = var[4]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(difcar, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(hyd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hyd2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(newcar, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(newhyd, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), difcar, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hyd1, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hyd2, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][0], p['HZ'][0], i), cal(p['0_H'][1], p['HZ'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', difcar), cal(p['CT_CT'][0], p['CT_CZ'][0], i), cal(p['CT_CT'][1], p['CT_CZ'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(difcar, hyd1), cal(p['CT_HC'][0], p['H_sCZ'][0], i), cal(p['CT_HC'][1], p['H_sCZ'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(difcar, hyd2), cal(p['CT_HC'][0], p['H_sCZ'][0], i), cal(p['CT_HC'][1], p['H_sCZ'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(difcar, newcar), cal(p['CT_mH'][0], p['CZ_CZ'][0], i), cal(p['CT_mH'][1], p['CZ_CZ'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, newhyd), cal(p['HC_mH'][0], p['CZ_HZ'][0], i), cal(p['HC_mH'][1], p['CZ_HZ'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', difcar), cal(p['C_C_C'][0], p['N_CT_C'][0], i), cal(p['C_C_C'][1], p['N_CT_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', difcar), cal(p['C_C_H'][0], p['H_C_CZ'][0], i), cal(p['C_C_H'][1], p['H_C_CZ'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', difcar, newcar), cal(p['C_C_H'][0], p['C_CZ_CZ'][0], i), cal(p['C_C_H'][1], p['C_CZ_CZ'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', difcar, hyd1), cal(p['C_C_H'][0], p['C_CZ_CZ'][0], i), cal(p['C_C_H'][1], p['C_CZ_CZ'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', difcar, hyd2), cal(p['C_C_H'][0], p['C_CZ_CZ'][0], i), cal(p['C_C_H'][1], p['C_CZ_CZ'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, difcar, hyd1), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, difcar, hyd2), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hyd1, difcar, hyd2), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hyd2, difcar, hyd2), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, newcar, difcar), cal(p['Dritt'][0], p['CZ_CZ_HZ'][0], i), cal(p['Dritt'][1], p['CZ_CZ_HZ'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hyd1, difcar, newcar, newhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hyd2, difcar, newcar, newhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', difcar, newcar, newhyd), cal(p['0_3'][0], p['0_3'][0], i), cal(p['0_3'][1], p['0_3'][1], i), cal(p['0_3'][2], p['0_3'][2], i), cal(p['0_3'][3], p['0_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', difcar, newcar), cal(p['X_CT_CZ_X'][0], p['X_CT_CZ_X'][0], i), cal(p['X_CT_CZ_X'][1], p['X_CT_CZ_X'][1], i), cal(p['X_CT_CZ_X'][2], p['X_CT_CZ_X'][2], i), cal(p['X_CT_CZ_X'][3], p['X_CT_CZ_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', difcar, hyd1), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', difcar, hyd2), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', difcar, newcar), cal(p['X_CT_CZ_X'][0], p['X_CT_CZ_X'][0], i), cal(p['X_CT_CZ_X'][1], p['X_CT_CZ_X'][1], i), cal(p['X_CT_CZ_X'][2], p['X_CT_CZ_X'][2], i), cal(p['X_CT_CZ_X'][3], p['X_CT_CZ_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', difcar, hyd1), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', difcar, hyd2), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), difcar, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hyd1, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hyd2, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][2], p['HZ'][2], i), cal(p['0_H'][3], p['HZ'][3], i))
