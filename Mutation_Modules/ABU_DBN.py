# ABU to DBN Mutation

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
        with open('Param_files/AminoAcid/DBN.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HG1 :{}@NE 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HG1 :{}@NZ 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@CG :{}@NZ 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vxi), bc['HG1']-(bc['HG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
		change(parm, 'charge', ':{}@ND'.format(vxi), fc['ND']/10*i*i/10).execute()
		change(parm, 'charge', ':{}@NE'.format(vxi), fc['NE']/10*i*i/10).execute()
		change(parm, 'charge', ':{}@NZ'.format(vxi), fc['NZ']/10*i*i/10).execute()
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
                        	pdb.write(atom.halfway_between1('ND', CG, HG1))
                        	pdb.write(atom.superimposed1('NE', HG1))
                        	pdb.write(atom.halfway_between2('NZ', CG, HG1))
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
	hydhyd = var[0]
	difhyd = var[1]
	azinit1 = var[2]
	azinit2 = var[3]
	azinit3 = var[4]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
        ctrl.write("%s=loadpdb Param_files/LibPDB/ABU-DBN.pdb\n"%vxi)
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
        ctrl.write('set %s.1.12 element "N"\n'%vxi)
        ctrl.write('set %s.1.13 element "N"\n'%vxi)
        ctrl.write('set %s.1.14 element "N"\n'%vxi)
        ctrl.write('set %s.1.15 element "C"\n'%vxi)
        ctrl.write('set %s.1.16 element "O"\n'%vxi)
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
        ctrl.write('set %s.1.12 name "ND"\n'%vxi)
        ctrl.write('set %s.1.13 name "NE"\n'%vxi)
        ctrl.write('set %s.1.14 name "NZ"\n'%vxi)
        ctrl.write('set %s.1.15 name "C"\n'%vxi)
        ctrl.write('set %s.1.16 name "O"\n'%vxi)
        ctrl.write('set %s.1.1 type "N"\n'%vxi)
        ctrl.write('set %s.1.2 type "H"\n'%vxi)
        ctrl.write('set %s.1.3 type "CT"\n'%vxi)
        ctrl.write('set %s.1.4 type "H1"\n'%vxi)
        ctrl.write('set %s.1.5 type "CT"\n'%vxi)
        ctrl.write('set %s.1.6 type "HC"\n'%vxi)
        ctrl.write('set %s.1.7 type "HC"\n'%vxi)
        ctrl.write('set %s.1.8 type "CT"\n'%vxi)
        ctrl.write('set %s.1.9 type "%s"\n'%(vxi, hydhyd))
        ctrl.write('set %s.1.10 type "%s"\n'%(vxi, difhyd))
        ctrl.write('set %s.1.11 type "%s"\n'%(vxi, difhyd))
        ctrl.write('set %s.1.12 type "%s"\n'%(vxi, azinit1))
        ctrl.write('set %s.1.13 type "%s"\n'%(vxi, azinit2))
        ctrl.write('set %s.1.14 type "%s"\n'%(vxi, azinit3))
        ctrl.write('set %s.1.15 type "C"\n'%vxi)
        ctrl.write('set %s.1.16 type "O"\n'%vxi)
        ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
        ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
        ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
        ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
        ctrl.write('bond %s.1.3 %s.1.15\n'%(vxi, vxi))
        ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
        ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
        ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
        ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
        ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
        ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
        ctrl.write('bond %s.1.8 %s.1.12\n'%(vxi, vxi))
        ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
        ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
        ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
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
	hydhyd = var[0]
	difhyd = var[1]
	azinit1 = var[2]
	azinit2 = var[3]
	azinit3 = var[4]
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(azinit1, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(azinit2, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(azinit3, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(hydhyd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(difhyd, 'H', 'sp3')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), difhyd, cal(p['HC'][0], p['HP'][0], i), cal(p['HC'][1], p['HP'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), azinit1, cal(p['NA'][0], p['NA'][0], i), cal(p['NA'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), azinit2, cal(p['0_N'][0], p['NA'][0], i), cal(p['0_N'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), azinit3, cal(p['0_N'][0], p['NA'][0], i), cal(p['0_N'][1], p['NA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', difhyd), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), cal(p['CT_HC'][0], p['H_sN'][0], i), cal(p['CT_HC'][1], p['H_sN'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', azinit1), cal(p['N_mHC'][0], p['C_NI'][0], i), cal(p['N_mHC'][1], p['C_NI'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(azinit1, azinit2), cal(p['N_mHC'][0], p['NI_ND'][0], i), cal(p['N_mHC'][1], p['NI_ND'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(azinit2, azinit3), cal(p['N_mHC'][0], p['ND_NE'][0], i), cal(p['N_mHC'][1], p['ND_NE'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd), cal(p['C_C_H'][0], p['C_C_NI'][0], i), cal(p['C_C_H'][1], p['C_C_NI'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', difhyd), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', azinit1), cal(p['C_C_N'][0], p['C_C_NI'][0], i), cal(p['C_C_N'][1], p['C_C_NI'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(difhyd, 'CT', azinit1), cal(p['H_C_H'][0], p['H_C_NI'][0], i), cal(p['H_C_H'][1], p['H_C_NI'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(difhyd, 'CT', hydhyd), cal(p['H_C_H'][0], p['H_C_NI'][0], i), cal(p['H_C_H'][1], p['H_C_NI'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(difhyd, 'CT', difhyd), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', azinit1), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', azinit1, azinit2), cal(p['Dritt'][0], p['C_NI_ND'][0], i), cal(p['Dritt'][1], p['C_NI_ND'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(azinit1, azinit2, azinit3), cal(p['Close'][0], p['NI_ND_NE'][0], i), cal(p['Close'][1], p['NI_ND_NE'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', azinit1, azinit2), cal(p['0_1'][0], p['X_C_NI_X_1'][0], i), cal(p['0_1'][1], p['X_C_NI_X_1'][1], i), cal(p['0_1'][2], p['X_C_NI_X_1'][2], i), cal(p['0_1'][3], p['X_C_NI_X_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', azinit1, azinit2), cal(p['0_2'][0], p['X_C_NI_X_2'][0], i), cal(p['0_2'][1], p['X_C_NI_X_2'][1], i), cal(p['0_2'][2], p['X_C_NI_X_2'][2], i), cal(p['0_2'][3], p['X_C_NI_X_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(difhyd, 'CT', azinit1, azinit2), cal(p['0_1'][0], p['X_C_NI_X_1'][0], i), cal(p['0_1'][1], p['X_C_NI_X_1'][1], i), cal(p['0_1'][2], p['X_C_NI_X_1'][2], i), cal(p['0_1'][3], p['X_C_NI_X_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(difhyd, 'CT', azinit1, azinit2), cal(p['0_2'][0], p['X_C_NI_X_2'][0], i), cal(p['0_2'][1], p['X_C_NI_X_2'][1], i), cal(p['0_2'][2], p['X_C_NI_X_2'][2], i), cal(p['0_2'][3], p['X_C_NI_X_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', azinit1, azinit2, azinit3), cal(p['0_3'][0], p['C_N_N_N'][0], i), cal(p['0_3'][1], p['C_N_N_N'][1], i), cal(p['0_3'][2], p['C_N_N_N'][2], i), cal(p['0_3'][3], p['C_N_N_N'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', azinit1, azinit2, azinit3), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), difhyd, cal(p['HC'][2], p['HP'][2], i), cal(p['HC'][3], p['HP'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), azinit1, cal(p['NA'][2], p['NA'][2], i), cal(p['NA'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), azinit2, cal(p['0_N'][2], p['NA'][2], i), cal(p['0_N'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), azinit3, cal(p['0_N'][2], p['NA'][2], i), cal(p['0_N'][3], p['NA'][3], i))
