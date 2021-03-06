# ALA to THR Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
	bc = {}
        with open('Param_files/AminoAcid/VAL.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/THR.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                changeLJPair(parm, ':{}@HB1 :{}@HG1 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HB :{}@HG11 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB'.format(vxi), fc['HB']/10*i).execute()
                change(parm, 'charge', ':{}@HB1'.format(vxi), bc['HB']-(bc['HB']/10)*i).execute()
                change(parm, 'charge', ':{}@CG1'.format(vxi), bc['CG1']-(bc['CG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG11'.format(vxi), bc['HG11']-(bc['HG11']/10)*i).execute()
                change(parm, 'charge', ':{}@HG12'.format(vxi), bc['HG12']-(bc['HG12']/10)*i).execute()
                change(parm, 'charge', ':{}@HG13'.format(vxi), bc['HG13']-(bc['HG13']/10)*i).execute()
                change(parm, 'charge', ':{}@CG2'.format(vxi), bc['CG2']+((fc['CG2']-bc['CG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG21'.format(vxi), bc['HG21']+((fc['HG21']-bc['HG21'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG22'.format(vxi), bc['HG22']+((fc['HG22']-bc['HG22'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG23'.format(vxi), bc['HG23']+((fc['HG23']-bc['HG23'])/10)*i).execute()
                change(parm, 'charge', ':{}@OG1'.format(vxi), (fc['OG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vxi), (fc['HG1']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CB = struct.residue_dict[aa].atom_dict['CB']
        HB = struct.residue_dict[aa].atom_dict['HB']
        CG1 = struct.residue_dict[aa].atom_dict['CG1']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HB' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HB1'))
                        	pdb.write(atom.superimposed1('HB', CG1))
			elif atom.get_name() == 'HG23' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('OG1', CB, HB))
                        	pdb.write(atom.superimposed1('HG1', HB))
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

def lib_make(ff, outputfile, vxi='VXI', metcar='dc', methyd='dh', hydhyd1='mh', alcoxy='ho', alchyd='hh', hydhyd2='sh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/VAL-THR.pdb\n"%vxi)
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
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "O"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB1"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG1"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG11"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG12"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG13"\n'%vxi)
	ctrl.write('set %s.1.12 name "CG2"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG21"\n'%vxi)
	ctrl.write('set %s.1.14 name "HG22"\n'%vxi)
	ctrl.write('set %s.1.15 name "HG23"\n'%vxi)
	ctrl.write('set %s.1.16 name "OG1"\n'%vxi)
	ctrl.write('set %s.1.17 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.18 name "C"\n'%vxi)
	ctrl.write('set %s.1.19 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, hydhyd1))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, hydhyd2))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.12 type "CT"\n'%vxi)
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "HC"\n'%vxi)
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, alcoxy))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, alchyd))
	ctrl.write('set %s.1.18 type "C"\n'%vxi)
	ctrl.write('set %s.1.19 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
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

def cal2(x, y, i):
        num = y+((x-y)/10)*i
        return num

def stock_add_to_all(metcar='dc', methyd='dh', hydhyd1='mh', alcoxy='ho', alchyd='hh', hydhyd2='sh', thrhyd='fh'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(metcar, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(methyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(alcoxy, 'O', 'sp3')
	Frcmod_creator.TYPE_insert(alchyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(thrhyd, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), alcoxy, cal2(p['OH'][0], p['0_O'][0], i), cal2(p['OH'][1], p['0_O'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), alchyd, cal2(p['HO'][0], p['0_H'][0], i), cal2(p['HO'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, cal2(p['0_H'][0], p['HC'][0], i), cal2(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), cal(p['CT_CT'][0], p['CT_mH'][0], i), cal(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', alcoxy), cal2(p['CT_OH'][0], p['OH_mH'][0], i), cal2(p['CT_OH'][1], p['OH_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd1), cal(p['HC_sC'][0], p['CT_HC'][0], i), cal(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd2), cal2(p['HC_sO'][0], p['CT_HC'][0], i), cal2(p['HC_sO'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), cal(p['CT_HC'][0], p['HC_mH'][0], i), cal(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(alcoxy, alchyd), cal2(p['OH_HO'][0], p['HO_mH'][0], i), cal2(p['OH_HO'][1], p['HO_mH'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', metcar), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', alcoxy), cal2(p['Close'][0], p['Close'][0], i), cal2(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', alcoxy), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', metcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', hydhyd2), cal2(p['H_C_H'][0], p['H_C_H'][0], i), cal2(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', alcoxy), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), cal(p['C_C_H'][0], p['Dritt'][0], i), cal(p['C_C_H'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', alcoxy, alchyd), cal2(p['C_O_H'][0], p['Dritt'][0], i), cal2(p['C_O_H'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(metcar, 'CT', alcoxy), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd1), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd2), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', metcar, methyd), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', alcoxy, alchyd), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', metcar, methyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', hydhyd1), cal2(p['X_C_O_X'][0], p['0_5'][0], i), cal2(p['X_C_O_X'][1], p['0_5'][1], i), cal2(p['X_C_O_X'][2], p['0_5'][2], i), cal2(p['X_C_O_X'][3], p['0_5'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar, 'CT', alcoxy, alchyd), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(methyd, metcar, 'CT', alcoxy), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', 'CT'), cal2(p['C_C_O_H_2'][0], p['0_3'][0], i), cal2(p['C_C_O_H_2'][1], p['0_3'][1], i), cal2(p['C_C_O_H_2'][2], p['0_3'][2], i), cal2(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', 'CT'), cal2(p['C_C_O_H_1'][0], p['0_2'][0], i), cal2(p['C_C_O_H_1'][1], p['0_2'][1], i), cal2(p['C_C_O_H_1'][2], p['0_2'][2], i), cal2(p['C_C_O_H_1'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', 'CT', 'H1'), cal2(p['C_C_O_H_2'][0], p['0_3'][0], i), cal2(p['C_C_O_H_2'][1], p['0_3'][1], i), cal2(p['C_C_O_H_2'][2], p['0_3'][2], i), cal2(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', 'CT', 'H1'), cal2(p['C_C_O_H_1'][0], p['0_2'][0], i), cal2(p['C_C_O_H_1'][1], p['0_2'][1], i), cal2(p['C_C_O_H_1'][2], p['0_2'][2], i), cal2(p['C_C_O_H_1'][3], p['0_2'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, cal(p['CT'][2], p['0_C'][2], i), cal(p['CT'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), alcoxy, cal2(p['OH'][2], p['0_O'][2], i), cal2(p['OH'][3], p['0_O'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), alchyd, cal2(p['HO'][2], p['0_H'][2], i), cal2(p['HO'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, cal2(p['0_H'][2], p['HC'][2], i), cal2(p['0_H'][3], p['HC'][3], i))
