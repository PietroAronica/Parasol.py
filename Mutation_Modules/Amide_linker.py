# Amide Stapler

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vx1='NX3', vx2='NX4'):
	bc = {}
        with open('Param_files/AminoAcid/ASN.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/AMX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HE22 :{}@HB1 0 0'.format(vx2, vx1)).execute()
                change(parm, 'charge', ':{}@N'.format(vx2), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx2), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx2), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx2), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx2), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx2), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx2), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx2), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@OD1'.format(vx2), bc['OD1']+((fc['OD1']-bc['OD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@ND2'.format(vx2), bc['ND2']+((fc['ND2']-bc['ND2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vx2), bc['HE21']+((fc['HE2']-bc['HE21'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE22'.format(vx2), bc['HE22']-(bc['HE22']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx2), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx2), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()
	bc = {}
        with open('Param_files/AminoAcid/ALA.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/ALX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB1'.format(vx1), bc['HB1']-(bc['HB1']/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, asn, ala, vx1='NX3', vx2='NX4'):
        struct.residue_dict[asn].set_resname(vx1)
        struct.residue_dict[ala].set_resname(vx2)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HE21' and res.get_resnumber() == asn:
                                pdb.write(atom.change_name('HE2'))
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

def all_make():
	for i in range(0,110,10):
		Frcmod_creator.make ('{}_{}.frcmod'.format(i, 100-i))

def cal(x, y, i):
        num = x+((y-x)/10)*i
        return num

def alc(x, y, i):
        num = x+((y-x)/10)*i*i/10
        return num

def lac(y, x, i):
        num = x+((y-x)/10)*i
        return num

def stock_add_to_all(dist, gonhyd='gh', dumhyd='dh'):
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobcar, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobhyd, cal(p['HC'][0], p['HA'][0], i), cal(p['HC'][1], p['HA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', dobcar), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, dobcar), alc(p['0_0'][0], p['CM_CM'][0], i), cal(dist, p['CM_CM'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, dobhyd), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, gonhyd), cal(p['CT_HC'][0], p['CA_HA2'][0], i), cal(p['CT_HC'][1], p['CA_HA2'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', dobcar), cal(p['C_C_C'][0], p['CT_CT_C'][0], i), cal(p['C_C_C'][1], p['CT_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', dobcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, dobhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, gonhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobhyd, dobcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['C_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, dobcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobcar, dobcar, gonhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobcar, dobcar, dobhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobcar, dobcar, 'CT'), cal(p['X_CM_CM_0'][0], p['C_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['C_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, dobhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_1'][0], p['0_1'][0], i), lac(p['C_C_C_CM_1'][1], p['0_1'][1], i), lac(p['C_C_C_CM_1'][2], p['0_1'][2], i), lac(p['C_C_C_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_2'][0], p['0_11'][0], i), lac(p['C_C_C_CM_2'][1], p['0_11'][1], i), lac(p['C_C_C_CM_2'][2], p['0_11'][2], i), lac(p['C_C_C_CM_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_3'][0], p['0_13'][0], i), lac(p['C_C_C_CM_3'][1], p['0_13'][1], i), lac(p['C_C_C_CM_3'][2], p['0_13'][2], i), lac(p['C_C_C_CM_3'][3], p['0_13'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, dobcar), lac(p['C_C_CM_CM_1'][0], p['0_1'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_1'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_1'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, dobcar), lac(p['C_C_CM_CM_2'][0], p['0_8'][0], i), lac(p['C_C_CM_CM_2'][1], p['0_8'][1], i), lac(p['C_C_CM_CM_2'][2], p['0_8'][2], i), lac(p['C_C_CM_CM_2'][3], p['0_8'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, dobcar), lac(p['C_C_CM_CM_3'][0], p['0_9'][0], i), lac(p['C_C_CM_CM_3'][1], p['0_9'][1], i), lac(p['C_C_CM_CM_3'][2], p['0_9'][2], i), lac(p['C_C_CM_CM_3'][3], p['0_9'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, dobcar), lac(p['H_C_CM_CM_1'][0], p['0_3'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_3'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_3'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, dobcar), lac(p['H_C_CM_CM_2'][0], p['0_14'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_14'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_14'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_14'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, dobhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, gonhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, dobhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, gonhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, dobcar, 'CT'), lac(p['C_CM_CM_C_1'][0], p['0_1'][0], i), lac(p['C_CM_CM_C_1'][1], p['0_1'][1], i), lac(p['C_CM_CM_C_1'][2], p['0_1'][2], i), lac(p['C_CM_CM_C_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, dobcar, 'CT'), lac(p['C_CM_CM_C_2'][0], p['0_11'][0], i), lac(p['C_CM_CM_C_2'][1], p['0_11'][1], i), lac(p['C_CM_CM_C_2'][2], p['0_11'][2], i), lac(p['C_CM_CM_C_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, dobcar, 'CT'), lac(p['C_CM_CM_C_3'][0], p['0_13'][0], i), lac(p['C_CM_CM_C_3'][1], p['0_13'][1], i), lac(p['C_CM_CM_C_3'][2], p['0_13'][2], i), lac(p['C_CM_CM_C_3'][3], p['0_13'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, dobcar, dobhyd), cal(p['0_15'][0], p['0_15'][0], i), cal(p['0_15'][1], p['0_15'][1], i), cal(p['0_15'][2], p['0_15'][2], i), cal(p['0_15'][3], p['0_15'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(dobhyd, dobcar, dobcar, dobhyd), cal(p['0_16'][0], p['0_16'][0], i), cal(p['0_16'][1], p['0_16'][1], i), cal(p['0_16'][2], p['0_16'][2], i), cal(p['0_16'][3], p['0_16'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, dobcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, dobcar, dobcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(dobhyd, dobcar, dobcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), dobcar, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), dobhyd, cal(p['HC'][2], p['HA'][2], i), cal(p['HC'][3], p['HA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
