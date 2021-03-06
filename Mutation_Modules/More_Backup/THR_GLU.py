# THR to GLU Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
	bc = {}
        with open('Param_files/AminoAcid/THR.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/GLU.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                changeLJPair(parm, ':{}@HB2 :{}@HG1 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HG21 :{}@OE1 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), fc['HB2']/10*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB']+((fc['HB3']-bc['HB'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG2']+((fc['CG']-bc['CG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG21'.format(vxi), bc['HG21']-(bc['HG21']/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG22']+((fc['HG2']-bc['HG22'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HG23']+((fc['HG3']-bc['HG23'])/10)*i).execute()
                change(parm, 'charge', ':{}@OG1'.format(vxi), bc['OG1']-(bc['OG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vxi), bc['HG1']-(bc['HG1']/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), (fc['CD']/10)*i*i/10).execute()
                change(parm, 'charge', ':{}@OE1'.format(vxi), (fc['OE1']/10)*i*i/10).execute()
                change(parm, 'charge', ':{}@OE2'.format(vxi), (fc['OE2']/10)*i*i/10).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		d = netCharge(parm).execute()
		change(parm, 'charge', ':PN', '{:.3f}'.format(-d)).execute() 
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CG2 = struct.residue_dict[aa].atom_dict['CG2']
        HG21 = struct.residue_dict[aa].atom_dict['HG21']
        OG1 = struct.residue_dict[aa].atom_dict['OG1']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HB' and res.get_resname() == vxi:
                        	pdb.write(atom.superimposed1('HB2', OG1))
                        	pdb.write(atom.change_name('HB3'))
			elif atom.get_name() == 'CG2' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CG'))
			elif atom.get_name() == 'HG22' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HG2'))
			elif atom.get_name() == 'HG23' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HG3'))
                        	pdb.write(atom.halfway_between('CD', CG2, HG21))
                        	pdb.write(atom.superimposed1('OE1', HG21))
                        	pdb.write(atom.superimposed2('OE2', HG21))
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

def lib_make(ff, outputfile, vxi='VXI', carcar='cc', caroxy='co', hydhyd1='yh', alcoxy='ho', alchyd='hh', hydhyd2='sh', thrhyd='fh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/THR-GLU.pdb\n"%vxi)
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
	ctrl.write('set %s.1.12 element "O"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "O"\n'%vxi)
	ctrl.write('set %s.1.16 element "O"\n'%vxi)
	ctrl.write('set %s.1.17 element "C"\n'%vxi)
	ctrl.write('set %s.1.18 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG21"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.12 name "OG1"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.14 name "CD"\n'%vxi)
	ctrl.write('set %s.1.15 name "OE1"\n'%vxi)
	ctrl.write('set %s.1.16 name "OE2"\n'%vxi)
	ctrl.write('set %s.1.17 name "C"\n'%vxi)
	ctrl.write('set %s.1.18 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, hydhyd2))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, thrhyd))
	ctrl.write('set %s.1.8 type "CT"\n'%vxi)
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, hydhyd1))
	ctrl.write('set %s.1.10 type "HC"\n'%vxi)
	ctrl.write('set %s.1.11 type "HC"\n'%vxi)
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, alcoxy))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, alchyd))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, carcar))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, caroxy))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, caroxy))
	ctrl.write('set %s.1.17 type "C"\n'%vxi)
	ctrl.write('set %s.1.18 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
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

def stock_add_to_all(vxi='VXI', carcar='cc', caroxy='co', hydhyd1='yh', alcoxy='ho', alchyd='hh', hydhyd2='sh', thrhyd='fh'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(alcoxy, 'O', 'sp3')
	Frcmod_creator.TYPE_insert(alchyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(carcar, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(caroxy, 'O', 'sp2')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), alcoxy, cal(p['OH'][0], p['0_O'][0], i), cal(p['OH'][1], p['0_O'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), alchyd, cal(p['HO'][0], p['0_H'][0], i), cal(p['HO'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), thrhyd, cal(p['H1'][0], p['HC'][0], i), cal(p['H1'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', alcoxy), cal(p['CT_OH'][0], p['OH_mH'][0], i), cal(p['CT_OH'][1], p['OH_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd2), cal(p['HC_sO'][0], p['CT_HC'][0], i), cal(p['HC_sO'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', thrhyd), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(alcoxy, alchyd), cal(p['OH_HO'][0], p['HO_mH'][0], i), cal(p['OH_HO'][1], p['HO_mH'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', alcoxy), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', alcoxy), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', alcoxy, alchyd), cal(p['C_O_H'][0], p['Dritt'][0], i), cal(p['C_O_H'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(thrhyd, 'CT', hydhyd2), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(thrhyd, 'CT', alcoxy), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', thrhyd), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd2), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', alcoxy, alchyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', thrhyd), cal(p['X_C_O_X'][0], p['0_5'][0], i), cal(p['X_C_O_X'][1], p['0_5'][1], i), cal(p['X_C_O_X'][2], p['0_5'][2], i), cal(p['X_C_O_X'][3], p['0_5'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', 'CT'), cal(p['C_C_O_H_2'][0], p['0_3'][0], i), cal(p['C_C_O_H_2'][1], p['0_3'][1], i), cal(p['C_C_O_H_2'][2], p['0_3'][2], i), cal(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', 'CT'), cal(p['C_C_O_H_1'][0], p['0_2'][0], i), cal(p['C_C_O_H_1'][1], p['0_2'][1], i), cal(p['C_C_O_H_1'][2], p['0_2'][2], i), cal(p['C_C_O_H_1'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', 'CT', 'H1'), cal(p['C_C_O_H_2'][0], p['0_3'][0], i), cal(p['C_C_O_H_2'][1], p['0_3'][1], i), cal(p['C_C_O_H_2'][2], p['0_3'][2], i), cal(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', 'CT', 'H1'), cal(p['C_C_O_H_1'][0], p['0_2'][0], i), cal(p['C_C_O_H_1'][1], p['0_2'][1], i), cal(p['C_C_O_H_1'][2], p['0_2'][2], i), cal(p['C_C_O_H_1'][3], p['0_2'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), alcoxy, cal(p['OH'][2], p['0_O'][2], i), cal(p['OH'][3], p['0_O'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), alchyd, cal(p['HO'][2], p['0_H'][2], i), cal(p['HO'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), thrhyd, cal(p['H1'][2], p['HC'][2], i), cal(p['H1'][3], p['HC'][3], i))

		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), carcar, cal2(p['C'][0], p['0_C'][0], i), cal2(p['C'][1], p['0_C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), caroxy, cal2(p['O2'][0], p['0_O'][0], i), cal2(p['O2'][1], p['0_O'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal2(p['0_H'][0], p['HC'][0], i), cal2(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', carcar), cal2(p['CT_C'][0], p['CT_mH'][0], i), cal2(p['CT_C'][1], p['CT_mH'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd1), cal2(p['HC_sC2'][0], p['CT_HC'][0], i), cal2(p['HC_sC2'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carcar, caroxy), cal2(p['C_O2'][0], p['O2_mH'][0], i), cal2(p['C_O2'][1], p['O2_mH'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', carcar, caroxy), cal2(p['C_C_O2'][0], p['Dritt'][0], i), cal2(p['C_C_O2'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(caroxy, carcar, caroxy), cal2(p['O2_C_O2'][0], p['Close'][0], i), cal2(p['O2_C_O2'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', carcar), cal2(p['CT_CT_C'][0], p['C_C_H'][0], i), cal2(p['CT_CT_C'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', carcar), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', hydhyd1), cal2(p['H_C_H'][0], p['H_C_H'][0], i), cal2(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', carcar), cal2(p['Close'][0], p['Close'][0], i), cal2(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd1), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carcar, caroxy), cal2(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal2(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal2(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal2(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carcar, caroxy), cal2(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal2(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal2(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal2(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', carcar, caroxy), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', caroxy, carcar, caroxy), cal2(p['Car_imp'][0], p['Imp_0'][0], i), cal2(p['Car_imp'][1], p['Imp_0'][1], i), cal2(p['Car_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), carcar, cal2(p['C'][2], p['0_C'][2], i), cal2(p['C'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), caroxy, cal2(p['O2'][2], p['0_O'][2], i), cal2(p['O2'][3], p['0_O'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal2(p['0_H'][2], p['HC'][2], i), cal2(p['0_H'][3], p['HC'][3], i))
