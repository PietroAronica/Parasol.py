# Remove GLY before ASP

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/ASP.param', 'r') as b:
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
		changeLJPair(parm, ':{}@CA :{}@H2 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HA :{}@H2 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@O :{}@H 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@O :{}@H1 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@O :{}@HA13 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N1'.format(vxi), (fc['N']-(fc['N']/10)*i)*(10-i)*(10-i)/100).execute()
                change(parm, 'charge', ':{}@H1'.format(vxi), (fc['H']-(fc['H']/10)*i)*(10-i)*(10-i)/100).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), fc['CA']+((bc['N']-fc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), fc['HA2']+((bc['H']-fc['HA2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA13'.format(vxi), (fc['HA3']-(fc['HA3']/10)*i)*(10-i)*(10-i)/100).execute()
                change(parm, 'charge', ':{}@C1'.format(vxi), (fc['C']-(fc['C']/10)*i)*(10-i)*(10-i)/100).execute()
                change(parm, 'charge', ':{}@O1'.format(vxi), (fc['O']-(fc['O']/10)*i)*(10-i)*(10-i)/100).execute()
                change(parm, 'charge', ':{}@N2'.format(vxi), (bc['N']-(bc['N']/10)*i)*(10-i)*(10-i)/100).execute()
                change(parm, 'charge', ':{}@H2'.format(vxi), (bc['H']-(bc['H']/10)*i)*(10-i)*(10-i)/100).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']).execute()
                change(parm, 'charge', ':{}@OD1'.format(vxi), bc['OD1']).execute()
                change(parm, 'charge', ':{}@OD2'.format(vxi), bc['OD2']).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']).execute()
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
			if atom.get_name() == 'N' and res.get_resnumber() == aa:
                        	pdb.write(atom.change_name('N1'))
			elif atom.get_name() == 'H' and res.get_resnumber() == aa:
                        	pdb.write(atom.change_name('H1'))
			elif atom.get_name() == 'CA' and res.get_resnumber() == aa:
                        	pdb.write(atom.change_name('N'))
			elif atom.get_name() == 'HA2' and res.get_resnumber() == aa:
                        	pdb.write(atom.change_name('H'))
			elif atom.get_name() == 'HA3' and res.get_resnumber() == aa:
                        	pdb.write(atom.change_name('HA13'))
			elif atom.get_name() == 'C' and res.get_resnumber() == aa:
                        	pdb.write(atom.change_name('C1'))
			elif atom.get_name() == 'O' and res.get_resnumber() == aa:
                        	pdb.write(atom.change_name('O1'))
			elif atom.get_resnum() == aa+1:
				atom.set_resname(vxi)
				atom.set_resnum(aa)
				if atom.get_name() == 'N':
                        		pdb.write(atom.change_name('N2'))
				elif atom.get_name() == 'H':
                        		pdb.write(atom.change_name('H2'))
				else:
	                        	pdb.write(atom.formatted())
			elif atom.get_resnum() > aa+1:
				atom.set_resnum(atom.get_resnum()-1)
                        	pdb.write(atom.formatted())
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
	var9 = sym + '9'
	var10 = sym + '0'
	var11 = sym + 'a'
	var12 = sym + 'b'
	var13 = sym + 'c'
	var14 = sym + 'd'
	var15 = sym + 'e'
	return var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, var15

def lib_make(ff, outputfile, vxi='VXI', var=variablemake()):
	oldnit = var[0]
	oldhyd = var[1]
	newnit = var[2]
	newhyd = var[3]
	newalp = var[4]
	newhl1 = var[5]
	newhl2 = var[6]
	newcar = var[7]
	newoxy = var[8]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/GLY_X_ASP.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "N"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "H"\n'%vxi)
	ctrl.write('set %s.1.6 element "C"\n'%vxi)
	ctrl.write('set %s.1.7 element "O"\n'%vxi)
	ctrl.write('set %s.1.8 element "N"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "O"\n'%vxi)
	ctrl.write('set %s.1.17 element "O"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N1"\n'%vxi)
	ctrl.write('set %s.1.2 name "H1"\n'%vxi)
	ctrl.write('set %s.1.3 name "N"\n'%vxi)
	ctrl.write('set %s.1.4 name "H"\n'%vxi)
	ctrl.write('set %s.1.5 name "HA13"\n'%vxi)
	ctrl.write('set %s.1.6 name "C1"\n'%vxi)
	ctrl.write('set %s.1.7 name "O1"\n'%vxi)
	ctrl.write('set %s.1.8 name "N2"\n'%vxi)
	ctrl.write('set %s.1.9 name "H2"\n'%vxi)
	ctrl.write('set %s.1.10 name "CA"\n'%vxi)
	ctrl.write('set %s.1.11 name "HA"\n'%vxi)
	ctrl.write('set %s.1.12 name "CB"\n'%vxi)
	ctrl.write('set %s.1.13 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.14 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.15 name "CG"\n'%vxi)
	ctrl.write('set %s.1.16 name "OD1"\n'%vxi)
	ctrl.write('set %s.1.17 name "OD2"\n'%vxi)
	ctrl.write('set %s.1.18 name "C"\n'%vxi)
	ctrl.write('set %s.1.19 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "%s"\n'%(vxi, newnit))
	ctrl.write('set %s.1.2 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.3 type "%s"\n'%(vxi, newalp))
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, newhl1))
	ctrl.write('set %s.1.5 type "%s"\n'%(vxi, newhl2))
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, newoxy))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, oldnit))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, oldhyd))
	ctrl.write('set %s.1.10 type "CX"\n'%vxi)
	ctrl.write('set %s.1.11 type "H1"\n'%vxi)
	ctrl.write('set %s.1.12 type "2C"\n'%vxi)
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "CO"\n'%vxi)
	ctrl.write('set %s.1.16 type "O2"\n'%vxi)
	ctrl.write('set %s.1.17 type "O2"\n'%vxi)
	ctrl.write('set %s.1.18 type "C"\n'%vxi)
	ctrl.write('set %s.1.19 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N1\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N1\n'%(vxi, vxi))
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
        num = float(y+((x-y)/10)*i)
        return float('{0:.3f}'.format(num))

def stock_add_to_all(var=variablemake()):
	oldnit = var[0]
	oldhyd = var[1]
	newnit = var[2]
	newhyd = var[3]
	newalp = var[4]
	newhl1 = var[5]
	newhl2 = var[6]
	newcar = var[7]
	newoxy = var[8]
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(oldnit, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(oldhyd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(newnit, 'N', 'sp3')
        Frcmod_creator.TYPE_insert(newhyd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(newalp, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(newhl1, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(newhl2, 'H', 'sp3')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newnit, lac(p['0_N'][0], p['NA'][0], i), lac(p['0_N'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newalp, lac(p['NA'][0], p['CT'][0], i), lac(p['NA'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhl1, lac(p['H'][0], p['H1'][0], i), lac(p['H'][1], p['H1'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhl2, lac(p['0_H'][0], p['H1'][0], i), lac(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, lac(p['0_C'][0], p['C'][0], i), lac(p['0_C'][1], p['C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newoxy, lac(p['0_O'][0], p['O2'][0], i), lac(p['0_O'][1], p['O2'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), oldnit, lac(p['0_N'][0], p['NA'][0], i), lac(p['0_N'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), oldhyd, lac(p['0_H'][0], p['H'][0], i), lac(p['0_H'][1], p['H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('C ', newnit), lac(p['Halve'][0], p['C_N'][0], i), lac(p['Halve'][1], p['C_N'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, newhyd), lac(p['Halve'][0], p['NA_H'][0], i), lac(p['Halve'][1], p['NA_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, newalp), lac(p['Halve'][0], p['CT_N2'][0], i), lac(p['Halve'][1], p['CT_N2'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newalp, newhl1), lac(p['NA_H'][0], p['CT_HC'][0], i), lac(p['NA_H'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newalp, newhl2), lac(p['NA_H'][0], p['CT_HC'][0], i), lac(p['NA_H'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newalp, newcar), lac(p['Third'][0], p['CT_C'][0], i), lac(p['Third'][1], p['CT_C'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, newoxy), lac(p['Third'][0], p['C_O'][0], i), lac(p['Third'][1], p['C_O'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, oldnit), lac(p['Third'][0], p['C_N'][0], i), lac(p['Third'][1], p['C_N'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(oldnit, oldhyd), lac(p['Third'][0], p['NA_H'][0], i), lac(p['Third'][1], p['NA_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(oldnit, 'CX'), lac(p['Third'][0], p['CT_N2'][0], i), lac(p['Third'][1], p['CT_N2'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('O ', 'C ', newnit), lac(p['O_C_N'][0], p['O_C_N'][0], i), lac(p['O_C_N'][1], p['O_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CX', 'C ', newnit), lac(p['C_C_N'][0], p['C_C_N'][0], i), lac(p['C_C_N'][1], p['C_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', newnit, newhyd), lac(p['Close'][0], p['F_CA_CA_HA'][0], i), lac(p['Close'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', newnit, newalp), lac(p['Dritt'][0], p['C_N_CT'][0], i), lac(p['Dritt'][1], p['C_N_CT'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, newnit, newalp), lac(p['Dritt'][0], p['CT_N_H'][0], i), lac(p['Dritt'][1], p['CT_N_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, newalp, newhl1), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, newalp, newhl2), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhl1, newalp, newcar), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhl2, newalp, newcar), lac(p['F_CA_CA_HA'][0], p['C_C_H'][0], i), lac(p['F_CA_CA_HA'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhl1, newalp, newhl2), lac(p['Close'][0], p['H_C_H'][0], i), lac(p['Close'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, newalp, newcar), lac(p['C_C_N'][0], p['N_CT_C'][0], i), lac(p['C_C_N'][1], p['N_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newalp, newcar, newoxy), lac(p['Close'][0], p['C_C_O'][0], i), lac(p['Close'][1], p['C_C_O'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newalp, newcar, oldnit), lac(p['Dritt'][0], p['C_C_N'][0], i), lac(p['Dritt'][1], p['C_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newoxy, newcar, oldnit), lac(p['Dritt'][0], p['O_C_N'][0], i), lac(p['Dritt'][1], p['O_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, oldnit, oldhyd), lac(p['Dritt'][0], p['F_CA_CA_HA'][0], i), lac(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, oldnit, 'CX'), lac(p['Dritt'][0], p['C_N_CT'][0], i), lac(p['Dritt'][1], p['C_N_CT'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(oldhyd, oldnit, 'CX'), lac(p['Close'][0], p['CT_N_H'][0], i), lac(p['Close'][1], p['CT_N_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(oldnit, 'CX', 'H1'), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(oldnit, 'CX', '2C'), lac(p['CT_CT_N'][0], p['CT_CT_N'][0], i), lac(p['CT_CT_N'][1], p['CT_CT_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(oldnit, 'CX', 'C '), lac(p['N_CT_C'][0], p['N_CT_C'][0], i), lac(p['N_CT_C'][1], p['N_CT_C'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CX', 'C ', newnit, newhyd), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('O ', 'C ', newnit, newhyd), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CX', 'C ', newnit, newalp), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('O ', 'C ', newnit, newalp), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, newalp, newhl1), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, newalp, newhl2), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, newalp, newcar), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, newalp, newcar, newoxy), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, newnit, newalp, newhl1), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, newnit, newalp, newhl2), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, newnit, newalp, newcar), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl1, newalp, newcar, newoxy), lac(p['0_10'][0], p['H_C_C_O_1'][0], i), lac(p['0_10'][1], p['H_C_C_O_1'][1], i), lac(p['0_10'][2], p['H_C_C_O_1'][2], i), lac(p['0_10'][3], p['H_C_C_O_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl1, newalp, newcar, newoxy), lac(p['0_8'][0], p['H_C_C_O_2'][0], i), lac(p['0_8'][1], p['H_C_C_O_2'][1], i), lac(p['0_8'][2], p['H_C_C_O_2'][2], i), lac(p['0_8'][3], p['H_C_C_O_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl1, newalp, newcar, newoxy), lac(p['0_9'][0], p['H_C_C_O_3'][0], i), lac(p['0_9'][1], p['H_C_C_O_3'][1], i), lac(p['0_9'][2], p['H_C_C_O_3'][2], i), lac(p['0_9'][3], p['H_C_C_O_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl2, newalp, newcar, newoxy), lac(p['0_10'][0], p['H_C_C_O_1'][0], i), lac(p['0_10'][1], p['H_C_C_O_1'][1], i), lac(p['0_10'][2], p['H_C_C_O_1'][2], i), lac(p['0_10'][3], p['H_C_C_O_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl2, newalp, newcar, newoxy), lac(p['0_8'][0], p['H_C_C_O_2'][0], i), lac(p['0_8'][1], p['H_C_C_O_2'][1], i), lac(p['0_8'][2], p['H_C_C_O_2'][2], i), lac(p['0_8'][3], p['H_C_C_O_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl2, newalp, newcar, newoxy), lac(p['0_9'][0], p['H_C_C_O_3'][0], i), lac(p['0_9'][1], p['H_C_C_O_3'][1], i), lac(p['0_9'][2], p['H_C_C_O_3'][2], i), lac(p['0_9'][3], p['H_C_C_O_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, newalp, newcar, oldnit), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl1, newalp, newcar, oldnit), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhl2, newalp, newcar, oldnit), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newalp, newcar, oldnit, oldhyd), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newalp, newcar, oldnit, 'CX'), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newoxy, newcar, oldnit, oldhyd), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newoxy, newcar, oldnit, 'CX'), lac(p['Ring_0'][0], p['X_C_N_X'][0], i), lac(p['Ring_0'][1], p['X_C_N_X'][1], i), lac(p['Ring_0'][2], p['X_C_N_X'][2], i), lac(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, oldnit, 'CX', 'H1'), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, oldnit, 'CX', '2C'), lac(p['C_C_C_N_1'][0], p['0_3'][0], i), lac(p['C_C_C_N_1'][1], p['0_3'][1], i), lac(p['C_C_C_N_1'][2], p['0_3'][2], i), lac(p['C_C_C_N_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, oldnit, 'CX', '2C'), lac(p['C_C_C_N_2'][0], p['0_8'][0], i), lac(p['C_C_C_N_2'][1], p['0_8'][1], i), lac(p['C_C_C_N_2'][2], p['0_8'][2], i), lac(p['C_C_C_N_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, oldnit, 'CX', '2C'), lac(p['C_C_C_N_3'][0], p['0_2'][0], i), lac(p['C_C_C_N_3'][1], p['0_2'][1], i), lac(p['C_C_C_N_3'][2], p['0_2'][2], i), lac(p['C_C_C_N_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, oldnit, 'CX', '2C'), lac(p['C_C_C_N_4'][0], p['0_7'][0], i), lac(p['C_C_C_N_4'][1], p['0_7'][1], i), lac(p['C_C_C_N_4'][2], p['0_7'][2], i), lac(p['C_C_C_N_4'][3], p['0_7'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, oldnit, 'CX', 'C '), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(oldhyd, oldnit, 'CX', 'H1'), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(oldhyd, oldnit, 'CX', '2C'), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(oldhyd, oldnit, 'CX', 'C '), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', oldnit, oldhyd), lac(p['Ami_imp'][0], p['Imp_0'][0], i), lac(p['Ami_imp'][1], p['Imp_0'][1], i), lac(p['Ami_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', newcar, newoxy), lac(p['Car_imp'][0], p['Imp_0'][0], i), lac(p['Car_imp'][1], p['Imp_0'][1], i), lac(p['Car_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', newnit, newhyd), lac(p['Ami_imp'][0], p['Imp_0'][0], i), lac(p['Ami_imp'][1], p['Imp_0'][1], i), lac(p['Ami_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newnit, lac(p['0_N'][2], p['NA'][2], i), lac(p['0_N'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newalp, lac(p['NA'][2], p['CT'][2], i), lac(p['NA'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhl1, lac(p['H'][2], p['H1'][2], i), lac(p['H'][3], p['H1'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhl2, lac(p['0_H'][2], p['H1'][2], i), lac(p['0_H'][3], p['H1'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, lac(p['0_C'][2], p['C'][2], i), lac(p['0_C'][3], p['C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newoxy, lac(p['0_O'][2], p['O2'][2], i), lac(p['0_O'][3], p['O2'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), oldnit, lac(p['0_N'][2], p['NA'][2], i), lac(p['0_N'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), oldhyd, lac(p['0_H'][2], p['H'][2], i), lac(p['0_H'][3], p['H'][3], i))
