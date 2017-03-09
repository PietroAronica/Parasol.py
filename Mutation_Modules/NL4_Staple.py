# NL4 Stapler

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='NLX'):
	bc = {}
        with open('Param_files/AminoAcid/NLE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/NL4.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HA :{}@HB21 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HE2 :{}@HE2 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']-(bc['HA']/10)*i).execute()
                change(parm, 'charge', ':{}@CB2'.format(vxi), (fc['CB2']/10)*i).execute()
                change(parm, 'charge', ':{}@HB21'.format(vxi), (fc['HB21']/10)*i).execute()
                change(parm, 'charge', ':{}@HB22'.format(vxi), (fc['HB22']/10)*i).execute()
                change(parm, 'charge', ':{}@HB23'.format(vxi), (fc['HB23']/10)*i).execute()
                change(parm, 'charge', ':{}@CB1'.format(vxi), bc['CB']+((fc['CB1']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB12'.format(vxi), bc['HB2']+((fc['HB12']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB13'.format(vxi), bc['HB3']+((fc['HB13']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), bc['CD']+((fc['CD']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), bc['HD3']+((fc['HD3']-bc['HD3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE'.format(vxi), bc['CE']+((fc['CE']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE'.format(vxi), bc['HE1']+((fc['HE']-bc['HE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), bc['HE2']-(bc['HE2']/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), bc['HE3']-(bc['HE3']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='NLX'):
	aa2 = aa + 4
        struct.residue_dict[aa].set_resname(vxi)
        struct.residue_dict[aa2].set_resname(vxi)
        CA1 = struct.residue_dict[aa].atom_dict['CA']
        HA1 = struct.residue_dict[aa].atom_dict['HA']
        CA2 = struct.residue_dict[aa2].atom_dict['CA']
        HA2 = struct.residue_dict[aa2].atom_dict['HA']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HA' and res.get_resnumber() == aa:
                        	pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('CB2', CA1, HA1))
                                pdb.write(atom.superimposed1('HB21', HA1))
                                pdb.write(atom.superimposed2('HB22', HA1))
                                pdb.write(atom.superimposed3('HB23', HA1))
			elif atom.get_name() == 'HA' and res.get_resnumber() == aa2:
                        	pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('CB2', CA2, HA2))
                                pdb.write(atom.superimposed1('HB21', HA2))
                                pdb.write(atom.superimposed2('HB22', HA2))
                                pdb.write(atom.superimposed3('HB23', HA2))
			elif atom.get_name() == 'CB' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('CB1'))
			elif atom.get_name() == 'HB2' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('HB12'))
			elif atom.get_name() == 'HB3' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('HB13'))
			elif atom.get_name() == 'HE1' and res.get_resname() == vxi:
                                pdb.write(atom.change_name('HE'))
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

def lib_make(ff, outputfile, vxi='VXI', metcar='dc', methyd='dh', hydhyd='mh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ABU-ALA.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB1"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.9 name "CG"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.12 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.13 name "C"\n'%vxi)
	ctrl.write('set %s.1.14 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, hydhyd))
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "HC"\n'%vxi)
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.13 type "C"\n'%vxi)
	ctrl.write('set %s.1.14 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.12\n'%(vxi, vxi))
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

def lac(y, x, i):
        num = x+((y-x)/10)*i
        return num

def stock_add_to_all(metcar='cb', methyd='hb', hydhyd='h1', dobcar='cm', dobhyd='ha', gonhyd='hc'):
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][0], p['0_C'][0], i), lac(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][0], p['H1'][0], i), lac(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), lac(p['CT_CT'][0], p['CT_mH'][0], i), lac(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), lac(p['HC_sC'][0], p['CT_HC'][0], i), lac(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), lac(p['CT_HC'][0], p['HC_mH'][0], i), lac(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), lac(p['C_C_H'][0], p['Dritt'][0], i), lac(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), lac(p['H_C_H'][0], p['Close'][0], i), lac(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', metcar), lac(p['CT_CT_C'][0], p['C_C_C'][0], i), lac(p['CT_CT_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', metcar), lac(p['CT_CT_N'][0], p['C_C_C'][0], i), lac(p['CT_CT_N'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), lac(p['C_C_C'][0], p['C_C_C'][0], i), lac(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', metcar), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_1'][0], p['0_3'][0], i), lac(p['C_C_N_C_1'][1], p['0_3'][1], i), lac(p['C_C_N_C_1'][2], p['0_3'][2], i), lac(p['C_C_N_C_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_2'][0], p['0_8'][0], i), lac(p['C_C_N_C_2'][1], p['0_8'][1], i), lac(p['C_C_N_C_2'][2], p['0_8'][2], i), lac(p['C_C_N_C_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_3'][0], p['0_2'][0], i), lac(p['C_C_N_C_3'][1], p['0_2'][1], i), lac(p['C_C_N_C_3'][2], p['0_2'][2], i), lac(p['C_C_N_C_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_4'][0], p['0_7'][0], i), lac(p['C_C_N_C_4'][1], p['0_7'][1], i), lac(p['C_C_N_C_4'][2], p['0_7'][2], i), lac(p['C_C_N_C_4'][3], p['0_7'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_1'][0], p['0_3'][0], i), lac(p['C_C_C_N_1'][1], p['0_3'][1], i), lac(p['C_C_C_N_1'][2], p['0_3'][2], i), lac(p['C_C_C_N_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_2'][0], p['0_8'][0], i), lac(p['C_C_C_N_2'][1], p['0_8'][1], i), lac(p['C_C_C_N_2'][2], p['0_8'][2], i), lac(p['C_C_C_N_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_3'][0], p['0_2'][0], i), lac(p['C_C_C_N_3'][1], p['0_2'][1], i), lac(p['C_C_C_N_3'][2], p['0_2'][2], i), lac(p['C_C_C_N_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_4'][0], p['0_7'][0], i), lac(p['C_C_C_N_4'][1], p['0_7'][1], i), lac(p['C_C_C_N_4'][2], p['0_7'][2], i), lac(p['C_C_C_N_4'][3], p['0_7'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), lac(p['C_C_C_H'][0], p['0_1'][0], i), lac(p['C_C_C_H'][1], p['0_1'][1], i), lac(p['C_C_C_H'][2], p['0_1'][2], i), lac(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', metcar, methyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'CT', metcar, methyd), lac(p['X_C_C_X'][0], p['0_4'][0], i), lac(p['X_C_C_X'][1], p['0_4'][1], i), lac(p['X_C_C_X'][2], p['0_4'][2], i), lac(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'CT', metcar, methyd), lac(p['X_C_C_X'][0], p['0_4'][0], i), lac(p['X_C_C_X'][1], p['0_4'][1], i), lac(p['X_C_C_X'][2], p['0_4'][2], i), lac(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][2], p['0_C'][2], i), lac(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][2], p['H1'][2], i), lac(p['0_H'][3], p['H1'][3], i))

                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobcar, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobhyd, cal(p['HC'][0], p['HA'][0], i), cal(p['HC'][1], p['HA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', dobcar), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, dobcar), cal(p['0_0'][0], p['CM_CM'][0], i), cal(p['0_0'][1], p['CM_CM'][1], i))
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
