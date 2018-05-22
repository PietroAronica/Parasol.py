# ALA to GLY Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/ALA.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/GLY.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HA2 :{}@HB1 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA2'.format(vxi), (fc['HA2']/10)*i).execute()
                change(parm, 'charge', ':{}@HA3'.format(vxi), bc['HA']+((fc['HA3']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']-(bc['CB']/10)*i).execute()
                change(parm, 'charge', ':{}@HB1'.format(vxi), bc['HB1']-(bc['HB1']/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']-(bc['HB2']/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']-(bc['HB3']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CB = struct.residue_dict[aa].atom_dict['CB']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'CA' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.superimposed1('HA2', CB))
			elif atom.get_name() == 'HA' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HA3'))
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

def lib_make(ff, outputfile, vxi='VXI', metcar='cd', methyd='dh', hydhyd='mh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ALA-GLY.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "H"\n'%vxi)
	ctrl.write('set %s.1.6 element "C"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA2"\n'%vxi)
	ctrl.write('set %s.1.5 name "HA3"\n'%vxi)
	ctrl.write('set %s.1.6 name "CB"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB1"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.9 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.10 name "C"\n'%vxi)
	ctrl.write('set %s.1.11 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, hydhyd))
	ctrl.write('set %s.1.5 type "H1"\n'%vxi)
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.10 type "C"\n'%vxi)
	ctrl.write('set %s.1.11 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.9\n'%(vxi, vxi))
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

def lac(x, y, i):
        num = y+((x-y)/10)*i
        return num

def stock_add_to_all(metcar='cd', methyd='dh', hydhyd='mh'):
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(metcar, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(methyd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(hydhyd, 'H', 'sp3')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, cal(p['CT'][0], p['0_C'][0], i), cal(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), cal(p['CT_CT'][0], p['CT_mH'][0], i), cal(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), cal(p['HC_sC'][0], p['CT_HC'][0], i), cal(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), cal(p['CT_HC'][0], p['HC_mH'][0], i), cal(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), cal(p['C_C_H'][0], p['Dritt'][0], i), cal(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', metcar), cal(p['CT_CT_C'][0], p['C_C_C'][0], i), cal(p['CT_CT_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', hydhyd), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', metcar), cal(p['CT_CT_N'][0], p['C_C_C'][0], i), cal(p['CT_CT_N'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', hydhyd), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', metcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', hydhyd), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', metcar), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), cal(p['C_C_N_C_1'][0], p['0_3'][0], i), cal(p['C_C_N_C_1'][1], p['0_3'][1], i), cal(p['C_C_N_C_1'][2], p['0_3'][2], i), cal(p['C_C_N_C_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), cal(p['C_C_N_C_2'][0], p['0_8'][0], i), cal(p['C_C_N_C_2'][1], p['0_8'][1], i), cal(p['C_C_N_C_2'][2], p['0_8'][2], i), cal(p['C_C_N_C_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), cal(p['C_C_N_C_3'][0], p['0_2'][0], i), cal(p['C_C_N_C_3'][1], p['0_2'][1], i), cal(p['C_C_N_C_3'][2], p['0_2'][2], i), cal(p['C_C_N_C_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), cal(p['C_C_N_C_4'][0], p['0_7'][0], i), cal(p['C_C_N_C_4'][1], p['0_7'][1], i), cal(p['C_C_N_C_4'][2], p['0_7'][2], i), cal(p['C_C_N_C_4'][3], p['0_7'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), cal(p['C_C_C_N_1'][0], p['0_3'][0], i), cal(p['C_C_C_N_1'][1], p['0_3'][1], i), cal(p['C_C_C_N_1'][2], p['0_3'][2], i), cal(p['C_C_C_N_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), cal(p['C_C_C_N_2'][0], p['0_8'][0], i), cal(p['C_C_C_N_2'][1], p['0_8'][1], i), cal(p['C_C_C_N_2'][2], p['0_8'][2], i), cal(p['C_C_C_N_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), cal(p['C_C_C_N_3'][0], p['0_2'][0], i), cal(p['C_C_C_N_3'][1], p['0_2'][1], i), cal(p['C_C_C_N_3'][2], p['0_2'][2], i), cal(p['C_C_C_N_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), cal(p['C_C_C_N_4'][0], p['0_7'][0], i), cal(p['C_C_C_N_4'][1], p['0_7'][1], i), cal(p['C_C_C_N_4'][2], p['0_7'][2], i), cal(p['C_C_C_N_4'][3], p['0_7'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', metcar, methyd), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', metcar, methyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'CT', metcar, methyd), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'CT', metcar, methyd), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, cal(p['CT'][2], p['0_C'][2], i), cal(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
