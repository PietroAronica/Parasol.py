# QUA to PRO Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/PRO.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/QUA.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG1'.format(vxi), bc['CG']+(bc['CG']/10)*i).execute()
                change(parm, 'charge', ':{}@HG12'.format(vxi), bc['HG2']+(bc['HG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HG13'.format(vxi), bc['HG3']+(bc['HG3']/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CD']+((fc['CG']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HD2']+((fc['HG2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HD3']+((fc['HG3']-bc['HD3'])/10)*i).execute()
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
                        if atom.get_name() == 'CG' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('CG1'))
                        elif atom.get_name() == 'HG2' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('HG12'))
                        elif atom.get_name() == 'HG3' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('HG13'))
                        elif atom.get_name() == 'CD' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('CG'))
                        elif atom.get_name() == 'HD2' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('HG2'))
                        elif atom.get_name() == 'HD3' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('HG3'))	
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

def lib_make(ff, outputfile, vxi='VXI', carcar='pc', procar='dc', prohyd='dh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/QUA-PRO.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "C"\n'%vxi)
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
	ctrl.write('set %s.1.14 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "CA"\n'%vxi)
	ctrl.write('set %s.1.3 name "HA"\n'%vxi)
	ctrl.write('set %s.1.4 name "CB"\n'%vxi)
	ctrl.write('set %s.1.5 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.7 name "CG1"\n'%vxi)
	ctrl.write('set %s.1.8 name "HG12"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG13"\n'%vxi)
	ctrl.write('set %s.1.10 name "CG"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.12 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.13 name "C"\n'%vxi)
	ctrl.write('set %s.1.14 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "CT"\n'%vxi)
	ctrl.write('set %s.1.3 type "H1"\n'%vxi)
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, carcar))
	ctrl.write('set %s.1.5 type "HC"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, procar))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, prohyd))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, prohyd))
	ctrl.write('set %s.1.10 type "CT"\n'%vxi)
	ctrl.write('set %s.1.11 type "H1"\n'%vxi)
	ctrl.write('set %s.1.12 type "H1"\n'%vxi)
	ctrl.write('set %s.1.13 type "C"\n'%vxi)
	ctrl.write('set %s.1.14 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.2 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.2 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.2 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
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

def stock_add_to_all(carcar='pc', procar='dc', prohyd='dh'):
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(carcar, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(procar, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(prohyd, 'H', 'sp3')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), carcar, lac(p['CT'][0], p['CT'][0], i), lac(p['CT'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), procar, lac(p['0_C'][0], p['CT'][0], i), lac(p['0_C'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), prohyd, lac(p['0_H'][0], p['HC'][0], i), lac(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carcar, procar), lac(p['CT_mCT'][0], p['CT_CT'][0], i), lac(p['CT_mCT'][1], p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', procar), lac(p['CT_mCT'][0], p['CT_CT'][0], i), lac(p['CT_mCT'][1], p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carcar, 'HC'), lac(p['CT_HC'][0], p['CT_HC'][0], i), lac(p['CT_HC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carcar, 'CT'), lac(p['CT_CT'][0], p['CT_CT'][0], i), lac(p['CT_CT'][1], p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(procar, prohyd), lac(p['HC_mCT'][0], p['CT_HC'][0], i), lac(p['HC_mCT'][1], p['CT_HC'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', procar, prohyd), lac(p['Dritt'][0], p['C_C_H'][0], i), lac(p['Dritt'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', procar, carcar), lac(p['Dritt'][0], p['C_C_C'][0], i), lac(p['Dritt'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carcar, procar, prohyd), lac(p['Close'][0], p['C_C_H'][0], i), lac(p['Close'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(prohyd, procar, prohyd), lac(p['Close'][0], p['H_C_H'][0], i), lac(p['Close'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', carcar, procar), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', carcar, 'HC'), lac(p['H_C_H'][0], p['H_C_H'][0], i), lac(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', carcar, 'HC'), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', carcar), lac(p['CT_CT_C'][0], p['CT_CT_C'][0], i), lac(p['CT_CT_C'][1], p['CT_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', carcar), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', carcar), lac(p['CT_CT_N'][0], p['CT_CT_N'][0], i), lac(p['CT_CT_N'][1], p['CT_CT_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', carcar, procar), lac(p['C_C_C'][0], p['C_C_C'][0], i), lac(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', procar), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', procar), lac(p['C_C_N'][0], p['C_C_N'][0], i), lac(p['C_C_N'][1], p['C_C_N'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N', 'CT', carcar, 'HC'), lac(p['0_4'][0], p['X_C_C_X'][0], i), lac(p['0_4'][1], p['X_C_C_X'][1], i), lac(p['0_4'][2], p['X_C_C_X'][2], i), lac(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C', 'CT', carcar, 'HC'), lac(p['0_4'][0], p['X_C_C_X'][0], i), lac(p['0_4'][1], p['X_C_C_X'][1], i), lac(p['0_4'][2], p['X_C_C_X'][2], i), lac(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', carcar, 'HC'), lac(p['0_1'][0], p['H_C_C_H'][0], i), lac(p['0_1'][1], p['H_C_C_H'][1], i), lac(p['0_1'][2], p['H_C_C_H'][2], i), lac(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N', 'CT', carcar, procar), lac(p['0_4'][0], p['X_C_C_X'][0], i), lac(p['0_4'][1], p['X_C_C_X'][1], i), lac(p['0_4'][2], p['X_C_C_X'][2], i), lac(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N', 'CT', procar, carcar), lac(p['0_4'][0], p['X_C_C_X'][0], i), lac(p['0_4'][1], p['X_C_C_X'][1], i), lac(p['0_4'][2], p['X_C_C_X'][2], i), lac(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N', 'CT', procar, prohyd), lac(p['0_1'][0], p['H_C_C_H'][0], i), lac(p['0_1'][1], p['H_C_C_H'][1], i), lac(p['0_1'][2], p['H_C_C_H'][2], i), lac(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C', 'CT', carcar, procar), lac(p['0_4'][0], p['X_C_C_X'][0], i), lac(p['0_4'][1], p['X_C_C_X'][1], i), lac(p['0_4'][2], p['X_C_C_X'][2], i), lac(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', carcar, procar), lac(p['0_1'][0], p['C_C_C_H'][0], i), lac(p['0_1'][1], p['C_C_C_H'][1], i), lac(p['0_1'][2], p['C_C_C_H'][2], i), lac(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', procar, prohyd), lac(p['0_1'][0], p['H_C_C_H'][0], i), lac(p['0_1'][1], p['H_C_C_H'][1], i), lac(p['0_1'][2], p['H_C_C_H'][2], i), lac(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', procar, carcar), lac(p['0_1'][0], p['C_C_C_H'][0], i), lac(p['0_1'][1], p['C_C_C_H'][1], i), lac(p['0_1'][2], p['C_C_C_H'][2], i), lac(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', carcar, procar, prohyd), lac(p['0_1'][0], p['H_C_C_H'][0], i), lac(p['0_1'][1], p['H_C_C_H'][1], i), lac(p['0_1'][2], p['H_C_C_H'][2], i), lac(p['0_1'][3], p['H_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', carcar, procar, prohyd), lac(p['0_1'][0], p['C_C_C_H'][0], i), lac(p['0_1'][1], p['C_C_C_H'][1], i), lac(p['0_1'][2], p['C_C_C_H'][2], i), lac(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', carcar, procar, 'CT'), lac(p['0_1'][0], p['C_C_C_H'][0], i), lac(p['0_1'][1], p['C_C_C_H'][1], i), lac(p['0_1'][2], p['C_C_C_H'][2], i), lac(p['0_1'][3], p['C_C_C_H'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', carcar, procar, 'CT'), lac(p['0_12'][0], p['C_C_C_C_1'][0], i), lac(p['0_12'][1], p['C_C_C_C_1'][1], i), lac(p['0_12'][2], p['C_C_C_C_1'][2], i), lac(p['0_12'][3], p['C_C_C_C_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', carcar, procar, 'CT'), lac(p['0_11'][0], p['C_C_C_C_2'][0], i), lac(p['0_11'][1], p['C_C_C_C_2'][1], i), lac(p['0_11'][2], p['C_C_C_C_2'][2], i), lac(p['0_11'][3], p['C_C_C_C_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', carcar, procar, 'CT'), lac(p['0_2'][0], p['C_C_C_C_3'][0], i), lac(p['0_2'][1], p['C_C_C_C_3'][1], i), lac(p['0_2'][2], p['C_C_C_C_3'][2], i), lac(p['0_2'][3], p['C_C_C_C_3'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), carcar, lac(p['CT'][2], p['CT'][2], i), lac(p['CT'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), procar, lac(p['0_C'][2], p['CT'][2], i), lac(p['0_C'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), prohyd, lac(p['0_H'][2], p['HC'][2], i), lac(p['0_H'][3], p['HC'][3], i))
