# Add Intercalated GLY to ALA

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
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
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']).execute()
                change(parm, 'charge', ':{}@HB1'.format(vxi), bc['HB1']).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']).execute()
                change(parm, 'charge', ':{}@N1'.format(vxi), (fc['N']/10)*i).execute()
                change(parm, 'charge', ':{}@H1'.format(vxi), (fc['H']/10)*i).execute()
                change(parm, 'charge', ':{}@CA1'.format(vxi), (fc['CA']/10)*i).execute()
                change(parm, 'charge', ':{}@HA12'.format(vxi), (fc['HA2']/10)*i).execute()
                change(parm, 'charge', ':{}@HA13'.format(vxi), (fc['HA3']/10)*i).execute()
                change(parm, 'charge', ':{}@C1'.format(vxi), (fc['C']/10)*i).execute()
                change(parm, 'charge', ':{}@O1'.format(vxi), (fc['O']/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        C = struct.residue_dict[aa].atom_dict['C']
        N = struct.residue_dict[aa+1].atom_dict['N']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'O' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.quarterone_between1('N1', C, N))
                        	pdb.write(atom.superimposed1('H1', C))
                        	pdb.write(atom.halfway_between('CA1', C, N))
                        	pdb.write(atom.quarterone_between2('HA12', C, N))
                        	pdb.write(atom.quarterone_between3('HA13', C, N))
                        	pdb.write(atom.quartertwo_between2('C1', C, N))
                        	pdb.write(atom.superimposed1('O1', N))
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

def lib_make(ff, outputfile, vxi='VXI', newnit='nn', newhyd='nh', newalp='c1', newhlp='h1', newcar='nc', newoxy='no'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ALA_X_GLY.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "O"\n'%vxi)
	ctrl.write('set %s.1.11 element "N"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB1"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.9 name "C"\n'%vxi)
	ctrl.write('set %s.1.10 name "O"\n'%vxi)
	ctrl.write('set %s.1.11 name "N1"\n'%vxi)
	ctrl.write('set %s.1.12 name "H1"\n'%vxi)
	ctrl.write('set %s.1.13 name "CA1"\n'%vxi)
	ctrl.write('set %s.1.14 name "HA12"\n'%vxi)
	ctrl.write('set %s.1.15 name "HA13"\n'%vxi)
	ctrl.write('set %s.1.16 name "C1"\n'%vxi)
	ctrl.write('set %s.1.17 name "O1"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "HC"\n'%vxi)
	ctrl.write('set %s.1.9 type "C"\n'%vxi)
	ctrl.write('set %s.1.10 type "O"\n'%vxi)
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, newnit))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, newalp))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, newhlp))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, newhlp))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, newoxy))
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C1\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C1\n'%(vxi, vxi))
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

def stock_add_to_all(newnit='nn', newhyd='nh', newalp='c1', newhlp='h1', newcar='nc', newoxy='no'):
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(newnit, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(newhyd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(newalp, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(newhlp, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(newcar, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(newoxy, 'O', 'sp3')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newnit, cal(p['0_N'][0], p['NA'][0], i), cal(p['0_N'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][0], p['H'][0], i), cal(p['0_H'][1], p['H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newalp, cal(p['0_C'][0], p['CT'][0], i), cal(p['0_C'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhlp, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['0_C'][0], p['C'][0], i), cal(p['0_C'][1], p['C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newoxy, cal(p['0_O'][0], p['O2'][0], i), cal(p['0_O'][1], p['O2'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('C ', newnit), cal(p['C_N_q'][0], p['C_N'][0], i), cal(p['C_N_q'][1], p['C_N'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, newhyd), cal(p['C_N_q'][0], p['NA_H'][0], i), cal(p['C_N_q'][1], p['NA_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, newalp), cal(p['C_N_q'][0], p['CT_N2'][0], i), cal(p['C_N_q'][1], p['CT_N2'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newalp, newhlp), cal(p['C_N_q'][0], p['CT_HC'][0], i), cal(p['C_N_q'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newalp, newcar), cal(p['C_N_q'][0], p['CT_C'][0], i), cal(p['C_N_q'][1], p['CT_C'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, newoxy), cal(p['C_N_q'][0], p['C_O'][0], i), cal(p['C_N_q'][1], p['C_O'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'N '), cal(p['C_N_q'][0], p['C_N'][0], i), cal(p['C_N_q'][1], p['C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('O ', 'C ', newnit), cal(p['O_C_N'][0], p['O_C_N'][0], i), cal(p['O_C_N'][1], p['O_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'C ', newnit), cal(p['C_C_N'][0], p['C_C_N'][0], i), cal(p['C_C_N'][1], p['C_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', newnit, newhyd), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', newnit, newalp), cal(p['Dritt'][0], p['C_N_CT'][0], i), cal(p['Dritt'][1], p['C_N_CT'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, newnit, newalp), cal(p['Close'][0], p['CT_N_H'][0], i), cal(p['Close'][1], p['CT_N_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, newalp, newhlp), cal(p['Close'][0], p['C_C_H'][0], i), cal(p['Close'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, newalp, newcar), cal(p['Dritt'][0], p['N_CT_C'][0], i), cal(p['Dritt'][1], p['N_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newalp, newcar, newoxy), cal(p['Dritt'][0], p['C_C_O'][0], i), cal(p['Dritt'][1], p['C_C_O'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhlp, newalp, newcar), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhlp, newalp, newhlp), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Dritt'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newalp, newcar, 'N '), cal(p['Dritt'][0], p['C_C_N'][0], i), cal(p['Dritt'][1], p['C_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newoxy, newcar, 'N '), cal(p['Close'][0], p['O_C_N'][0], i), cal(p['Close'][1], p['O_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'N ', 'H '), cal(p['F_CA_CA_HA'][0], p['F_CA_CA_HA'][0], i), cal(p['F_CA_CA_HA'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'N ', 'CX'), cal(p['C_N_CT'][0], p['C_N_CT'][0], i), cal(p['C_N_CT'][1], p['C_N_CT'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'C ', newnit, newhyd), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'C ', newnit, newalp), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('O ', 'C ', newnit, newhyd), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('O ', 'C ', newnit, newalp), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, newalp, newhlp), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, newalp, newcar), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, newnit, newalp, newhlp), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, newnit, newalp, newcar), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhlp, newalp, newcar, newoxy), cal(p['0_10'][0], p['H_C_C_O_1'][0], i), cal(p['0_10'][1], p['H_C_C_O_1'][1], i), cal(p['0_10'][2], p['H_C_C_O_1'][2], i), cal(p['0_10'][3], p['H_C_C_O_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhlp, newalp, newcar, newoxy), cal(p['0_8'][0], p['H_C_C_O_2'][0], i), cal(p['0_8'][1], p['H_C_C_O_2'][1], i), cal(p['0_8'][2], p['H_C_C_O_2'][2], i), cal(p['0_8'][3], p['H_C_C_O_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhlp, newalp, newcar, newoxy), cal(p['0_9'][0], p['H_C_C_O_3'][0], i), cal(p['0_9'][1], p['H_C_C_O_3'][1], i), cal(p['0_9'][2], p['H_C_C_O_3'][2], i), cal(p['0_9'][3], p['H_C_C_O_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, newalp, newcar, newoxy), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, newalp, newcar, 'N '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhlp, newalp, newcar, 'N '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newalp, newcar, 'N ', 'H '), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newoxy, newcar, 'N ', 'H '), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newalp, newcar, 'N ', 'CX'), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newoxy, newcar, 'N ', 'CX'), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', newcar, newoxy), lac(p['Car_imp'][0], p['Imp_0'][0], i), lac(p['Car_imp'][1], p['Imp_0'][1], i), lac(p['Car_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', newnit, newhyd), lac(p['Ami_imp'][0], p['Imp_0'][0], i), lac(p['Ami_imp'][1], p['Imp_0'][1], i), lac(p['Ami_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newnit, cal(p['0_N'][2], p['NA'][2], i), cal(p['0_N'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['0_H'][2], p['H'][2], i), cal(p['0_H'][3], p['H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newalp, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhlp, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['0_C'][2], p['C'][2], i), cal(p['0_C'][3], p['C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newoxy, cal(p['0_O'][2], p['O2'][2], i), cal(p['0_O'][3], p['O2'][3], i))
