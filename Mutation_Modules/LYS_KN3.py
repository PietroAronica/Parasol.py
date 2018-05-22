# LYS to KN3 Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/LYS.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/KN3.param', 'r') as b:
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
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), bc['CD']+((fc['CD']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vxi), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vxi), bc['HD3']+((fc['HD3']-bc['HD3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE'.format(vxi), bc['CE']+((fc['CE']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), bc['HE2']+((fc['HE2']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), bc['HE3']+((fc['HE3']-bc['HE3'])/10)*i).execute()
                change(parm, 'charge', ':{}@NZ'.format(vxi), bc['NZ']-(bc['NZ']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ1'.format(vxi), bc['HZ1']-(bc['HZ1']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vxi), bc['HZ2']-(bc['HZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@HZ3'.format(vxi), bc['HZ3']-(bc['HZ3']/10)*i).execute()
                change(parm, 'charge', ':{}@NH'.format(vxi), (fc['NH']/10)*i).execute()
                change(parm, 'charge', ':{}@NT'.format(vxi), (fc['NT']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
                d = netCharge(parm).execute()
                change(parm, 'charge', ':PN', '{:.3f}'.format(-d)).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        NZ = struct.residue_dict[aa].atom_dict['NZ']
        HZ1 = struct.residue_dict[aa].atom_dict['HZ1']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HZ3' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.superimposed1('NH', HZ1))
                        	pdb.write(atom.superimposed1('NT', NZ))
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

def lib_make(ff, outputfile, vxi='VXI', hyd1='h1', hyd2='h2', azinit1='n1', azinit2='n2', azinit3='n3'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
        ctrl.write("%s=loadpdb Param_files/LibPDB/LYS-KN3.pdb\n"%vxi)
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
        ctrl.write('set %s.1.11 element "C"\n'%vxi)
        ctrl.write('set %s.1.12 element "H"\n'%vxi)
        ctrl.write('set %s.1.13 element "H"\n'%vxi)
        ctrl.write('set %s.1.14 element "C"\n'%vxi)
        ctrl.write('set %s.1.15 element "H"\n'%vxi)
        ctrl.write('set %s.1.16 element "H"\n'%vxi)
        ctrl.write('set %s.1.17 element "N"\n'%vxi)
        ctrl.write('set %s.1.18 element "H"\n'%vxi)
        ctrl.write('set %s.1.19 element "H"\n'%vxi)
        ctrl.write('set %s.1.20 element "H"\n'%vxi)
        ctrl.write('set %s.1.21 element "N"\n'%vxi)
        ctrl.write('set %s.1.22 element "N"\n'%vxi)
        ctrl.write('set %s.1.23 element "C"\n'%vxi)
        ctrl.write('set %s.1.24 element "O"\n'%vxi)
        ctrl.write('set %s.1.1 name "N"\n'%vxi)
        ctrl.write('set %s.1.2 name "H"\n'%vxi)
        ctrl.write('set %s.1.3 name "CA"\n'%vxi)
        ctrl.write('set %s.1.4 name "HA"\n'%vxi)
        ctrl.write('set %s.1.5 name "CB"\n'%vxi)
        ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
        ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
        ctrl.write('set %s.1.8 name "CG"\n'%vxi)
        ctrl.write('set %s.1.9 name "HG2"\n'%vxi)
        ctrl.write('set %s.1.10 name "HG3"\n'%vxi)
        ctrl.write('set %s.1.11 name "CD"\n'%vxi)
        ctrl.write('set %s.1.12 name "HD2"\n'%vxi)
        ctrl.write('set %s.1.13 name "HD3"\n'%vxi)
        ctrl.write('set %s.1.14 name "CE"\n'%vxi)
        ctrl.write('set %s.1.15 name "HE2"\n'%vxi)
        ctrl.write('set %s.1.16 name "HE3"\n'%vxi)
        ctrl.write('set %s.1.17 name "NZ"\n'%vxi)
        ctrl.write('set %s.1.18 name "HZ1"\n'%vxi)
        ctrl.write('set %s.1.19 name "HZ2"\n'%vxi)
        ctrl.write('set %s.1.20 name "HZ3"\n'%vxi)
        ctrl.write('set %s.1.21 name "NH"\n'%vxi)
        ctrl.write('set %s.1.22 name "NT"\n'%vxi)
        ctrl.write('set %s.1.23 name "C"\n'%vxi)
        ctrl.write('set %s.1.24 name "O"\n'%vxi)
        ctrl.write('set %s.1.1 type "N"\n'%vxi)
        ctrl.write('set %s.1.2 type "H"\n'%vxi)
        ctrl.write('set %s.1.3 type "CT"\n'%vxi)
        ctrl.write('set %s.1.4 type "H1"\n'%vxi)
        ctrl.write('set %s.1.5 type "CT"\n'%vxi)
        ctrl.write('set %s.1.6 type "HC"\n'%vxi)
        ctrl.write('set %s.1.7 type "HC"\n'%vxi)
        ctrl.write('set %s.1.8 type "CT"\n'%vxi)
        ctrl.write('set %s.1.9 type "HC"\n'%vxi)
        ctrl.write('set %s.1.10 type "HC"\n'%vxi)
        ctrl.write('set %s.1.11 type "CT"\n'%vxi)
        ctrl.write('set %s.1.12 type "HC"\n'%vxi)
        ctrl.write('set %s.1.13 type "HC"\n'%vxi)
        ctrl.write('set %s.1.14 type "CT"\n'%vxi)
        ctrl.write('set %s.1.15 type "HP"\n'%vxi)
        ctrl.write('set %s.1.16 type "HP"\n'%vxi)
        ctrl.write('set %s.1.17 type "%s"\n'%(vxi, azinit1))
        ctrl.write('set %s.1.18 type "%s"\n'%(vxi, hyd1))
        ctrl.write('set %s.1.19 type "%s"\n'%(vxi, hyd2))
        ctrl.write('set %s.1.20 type "%s"\n'%(vxi, hyd2))
        ctrl.write('set %s.1.21 type "%s"\n'%(vxi, azinit2))
        ctrl.write('set %s.1.22 type "%s"\n'%(vxi, azinit3))
        ctrl.write('set %s.1.23 type "C"\n'%vxi)
        ctrl.write('set %s.1.24 type "O"\n'%vxi)
        ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
        ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
        ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
        ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
        ctrl.write('bond %s.1.3 %s.1.23\n'%(vxi, vxi))
        ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
        ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
        ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
        ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
        ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
        ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
        ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
        ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
        ctrl.write('bond %s.1.11 %s.1.14\n'%(vxi, vxi))
        ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
        ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
        ctrl.write('bond %s.1.14 %s.1.17\n'%(vxi, vxi))
        ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
        ctrl.write('bond %s.1.17 %s.1.19\n'%(vxi, vxi))
        ctrl.write('bond %s.1.17 %s.1.20\n'%(vxi, vxi))
        ctrl.write('bond %s.1.17 %s.1.21\n'%(vxi, vxi))
        ctrl.write('bond %s.1.21 %s.1.22\n'%(vxi, vxi))
        ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
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

def stock_add_to_all(hyd1='h1', hyd2='h2', azinit1='n1', azinit2='n2', azinit3='n3'):
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(azinit1, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(azinit2, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(azinit3, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(hyd1, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(hyd2, 'H', 'sp3')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), azinit1, cal(p['NA'][0], p['NA'][0], i), cal(p['NA'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), azinit2, cal(p['0_N'][0], p['NA'][0], i), cal(p['0_N'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), azinit3, cal(p['0_N'][0], p['NA'][0], i), cal(p['0_N'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hyd1, cal(p['H'][0], p['0_H'][0], i), cal(p['H'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hyd2, cal(p['H'][0], p['0_H'][0], i), cal(p['H'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', azinit1), cal(p['C_N3'][0], p['C_NI'][0], i), cal(p['C_N3'][1], p['C_NI'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(azinit1, hyd1), cal(p['NA_H'][0], p['NI_sH'][0], i), cal(p['NA_H'][1], p['NI_sH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(azinit1, hyd2), cal(p['NA_H'][0], p['NI_sH'][0], i), cal(p['NA_H'][1], p['NI_sH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(azinit1, azinit2), cal(p['ND_sN'][0], p['NI_ND'][0], i), cal(p['ND_sN'][1], p['NI_ND'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(azinit2, azinit3), cal(p['NE_sN'][0], p['ND_NE'][0], i), cal(p['NE_sN'][1], p['ND_NE'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', azinit1), cal(p['C_C_N'][0], p['C_C_NI'][0], i), cal(p['C_C_N'][1], p['C_C_NI'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HP', 'CT', azinit1), cal(p['C_C_H'][0], p['H_C_NI'][0], i), cal(p['C_C_H'][1], p['H_C_NI'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', azinit1, hyd1), cal(p['C_C_H'][0], p['N_sH'][0], i), cal(p['C_C_H'][1], p['N_sH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', azinit1, hyd2), cal(p['C_C_H'][0], p['N_sH'][0], i), cal(p['C_C_H'][1], p['N_sH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', azinit1, azinit2), cal(p['N_sH'][0], p['C_NI_ND'][0], i), cal(p['N_sH'][1], p['C_NI_ND'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(azinit1, azinit2, azinit3), cal(p['Close'][0], p['NI_ND_NE'][0], i), cal(p['Close'][1], p['NI_ND_NE'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hyd1, azinit1, azinit2), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hyd2, azinit1, azinit2), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hyd1, azinit1, hyd2), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hyd2, azinit1, hyd2), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', azinit1, hyd1), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', azinit1, hyd2), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HP', 'CT', azinit1, hyd1), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HP', 'CT', azinit1, hyd2), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', azinit1, azinit2), cal(p['0_1'][0], p['X_C_NI_X_1'][0], i), cal(p['0_1'][1], p['X_C_NI_X_1'][1], i), cal(p['0_1'][2], p['X_C_NI_X_1'][2], i), cal(p['0_1'][3], p['X_C_NI_X_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', azinit1, azinit2), cal(p['0_2'][0], p['X_C_NI_X_2'][0], i), cal(p['0_2'][1], p['X_C_NI_X_2'][1], i), cal(p['0_2'][2], p['X_C_NI_X_2'][2], i), cal(p['0_2'][3], p['X_C_NI_X_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HP', 'CT', azinit1, azinit2), cal(p['0_1'][0], p['X_C_NI_X_1'][0], i), cal(p['0_1'][1], p['X_C_NI_X_1'][1], i), cal(p['0_1'][2], p['X_C_NI_X_1'][2], i), cal(p['0_1'][3], p['X_C_NI_X_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HP', 'CT', azinit1, azinit2), cal(p['0_2'][0], p['X_C_NI_X_2'][0], i), cal(p['0_2'][1], p['X_C_NI_X_2'][1], i), cal(p['0_2'][2], p['X_C_NI_X_2'][2], i), cal(p['0_2'][3], p['X_C_NI_X_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', azinit1, azinit2, azinit3), cal(p['0_3'][0], p['C_N_N_N'][0], i), cal(p['0_3'][1], p['C_N_N_N'][1], i), cal(p['0_3'][2], p['C_N_N_N'][2], i), cal(p['0_3'][3], p['C_N_N_N'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hyd1, azinit1, azinit2, azinit3), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hyd2, azinit1, azinit2, azinit3), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), azinit1, cal(p['NA'][2], p['NA'][2], i), cal(p['NA'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), azinit2, cal(p['0_N'][2], p['NA'][2], i), cal(p['0_N'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), azinit3, cal(p['0_N'][2], p['NA'][2], i), cal(p['0_N'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hyd1, cal(p['H'][2], p['0_H'][2], i), cal(p['H'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hyd2, cal(p['H'][2], p['0_H'][2], i), cal(p['H'][3], p['0_H'][3], i))
