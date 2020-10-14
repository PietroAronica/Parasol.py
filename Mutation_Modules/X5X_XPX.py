# X5X to XPX Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/X5X.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/XPX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@H32'.format(vxi), ':{}@H22'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@H33'.format(vxi), ':{}@H22'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@C21'.format(vxi), ':{}@H22'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@C21'.format(vxi), ':{}@C32'.format(vxi), '0', '0').execute()
                change(parm, 'charge', ':{}@CA1'.format(vxi), bc['C1']+((fc['CA1']-bc['C1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA12'.format(vxi), bc['H12']+((fc['HA12']-bc['H12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA13'.format(vxi), bc['H13']+((fc['HA13']-bc['H13'])/10)*i).execute()
                change(parm, 'charge', ':{}@C1'.format(vxi), bc['C2']+((fc['C1']-bc['C2'])/10)*i).execute()
                change(parm, 'charge', ':{}@H12'.format(vxi), bc['H22']-(bc['H22']/10)*i).execute()
                change(parm, 'charge', ':{}@H13'.format(vxi), bc['H23']-(bc['H23']/10)*i).execute()
                change(parm, 'charge', ':{}@C21'.format(vxi), (fc['C21']/10)*i).execute()
                change(parm, 'charge', ':{}@H21'.format(vxi), (fc['H21']/10)*i).execute()
                change(parm, 'charge', ':{}@C31'.format(vxi), (fc['C31']/10)*i).execute()
                change(parm, 'charge', ':{}@H31'.format(vxi), (fc['H31']/10)*i).execute()
                change(parm, 'charge', ':{}@C4'.format(vxi), bc['C4']+((fc['C4']-bc['C4'])/10)*i).execute()
                change(parm, 'charge', ':{}@H42'.format(vxi), bc['H42']-(bc['H42']/10)*i).execute()
                change(parm, 'charge', ':{}@H43'.format(vxi), bc['H43']-(bc['H43']/10)*i).execute()
                change(parm, 'charge', ':{}@C32'.format(vxi), bc['C3']+((fc['C32']-bc['C3'])/10)*i).execute()
                change(parm, 'charge', ':{}@H32'.format(vxi), bc['H32']+((fc['H32']-bc['H32'])/10)*i).execute()
                change(parm, 'charge', ':{}@H33'.format(vxi), bc['H33']-(bc['H33']/10)*i).execute()
                change(parm, 'charge', ':{}@C22'.format(vxi), (fc['C22']/10)*i).execute()
                change(parm, 'charge', ':{}@H22'.format(vxi), (fc['H22']/10)*i).execute()
                change(parm, 'charge', ':{}@CA2'.format(vxi), bc['C5']+((fc['CA2']-bc['C5'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA22'.format(vxi), bc['H52']+((fc['HA22']-bc['H52'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA23'.format(vxi), bc['H53']+((fc['HA23']-bc['H53'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        C2 = struct.residue_dict[aa].atom_dict['C2']
        C3 = struct.residue_dict[aa].atom_dict['C3']
        C4 = struct.residue_dict[aa].atom_dict['C4']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'C1' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CA1'))
			elif atom.get_name() == 'H12' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HA12'))
			elif atom.get_name() == 'H13' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HA13'))
			elif atom.get_name() == 'C2' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('C1'))
			elif atom.get_name() == 'H22' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('H12'))
			elif atom.get_name() == 'H23' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('H13'))
                        	pdb.write(atom.superimposed1('C21', C3))
                        	pdb.write(atom.superimposed1('H21', C2))
                        	pdb.write(atom.halfway_between('C31', C3, C4))
                        	pdb.write(atom.superimposed1('H31', C4))
			elif atom.get_name() == 'C3' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('C32'))
			elif atom.get_name() == 'H33' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('C22', C2, C3))
                        	pdb.write(atom.superimposed2('H22', C3))
			elif atom.get_name() == 'C5' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CA2'))
			elif atom.get_name() == 'H52' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HA22'))
			elif atom.get_name() == 'H53' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HA23'))
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
	return var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14

def lib_make(ff, outputfile, vxi='VXI', var=variablemake()):
 	c1 = var[0]
 	h1 = var[1]
        c21 = var[2]
        h21 = var[3]
        c31 = var[4]
        h31 = var[5]
        c4 = var[6]
        h4 = var[7]
        c32 = var[8]
        h32 = var[9]
        h3 = var[10]
        c22 = var[11]
        h22 = var[12]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/X5X-XPX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "C"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "H"\n'%vxi)
	ctrl.write('set %s.1.4 element "C"\n'%vxi)
	ctrl.write('set %s.1.5 element "H"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "C"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "C"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.1 name "CA1"\n'%vxi)
	ctrl.write('set %s.1.2 name "HA12"\n'%vxi)
	ctrl.write('set %s.1.3 name "HA13"\n'%vxi)
	ctrl.write('set %s.1.4 name "C1"\n'%vxi)
	ctrl.write('set %s.1.5 name "H12"\n'%vxi)
	ctrl.write('set %s.1.6 name "H13"\n'%vxi)
	ctrl.write('set %s.1.7 name "C21"\n'%vxi)
	ctrl.write('set %s.1.8 name "H21"\n'%vxi)
	ctrl.write('set %s.1.9 name "C31"\n'%vxi)
	ctrl.write('set %s.1.10 name "H31"\n'%vxi)
	ctrl.write('set %s.1.11 name "C4"\n'%vxi)
	ctrl.write('set %s.1.12 name "H42"\n'%vxi)
	ctrl.write('set %s.1.13 name "H43"\n'%vxi)
	ctrl.write('set %s.1.14 name "C32"\n'%vxi)
	ctrl.write('set %s.1.15 name "H32"\n'%vxi)
	ctrl.write('set %s.1.16 name "H33"\n'%vxi)
	ctrl.write('set %s.1.17 name "C22"\n'%vxi)
	ctrl.write('set %s.1.18 name "H22"\n'%vxi)
	ctrl.write('set %s.1.19 name "CA2"\n'%vxi)
	ctrl.write('set %s.1.20 name "HA22"\n'%vxi)
	ctrl.write('set %s.1.21 name "HA23"\n'%vxi)
	ctrl.write('set %s.1.1 type "CT"\n'%vxi)
	ctrl.write('set %s.1.2 type "HC"\n'%vxi)
	ctrl.write('set %s.1.3 type "HC"\n'%vxi)
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, c1))
	ctrl.write('set %s.1.5 type "%s"\n'%(vxi, h1))
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, h1))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, c21))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, h21))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, c31))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, h31))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, c4))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, h4))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, h4))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, c32))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, h32))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, h3))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, c22))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, h22))
	ctrl.write('set %s.1.19 type "CT"\n'%vxi)
	ctrl.write('set %s.1.20 type "HC"\n'%vxi)
	ctrl.write('set %s.1.21 type "HC"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.17 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.21\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.CA1\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.CA2\n'%(vxi, vxi))
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
 	c1 = var[0]
 	h1 = var[1]
        c21 = var[2]
        h21 = var[3]
        c31 = var[4]
        h31 = var[5]
        c4 = var[6]
        h4 = var[7]
        c32 = var[8]
        h32 = var[9]
        h3 = var[10]
        c22 = var[11]
        h22 = var[12]
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(c1, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(h1, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(c21, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(h21, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(c31, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(h31, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(c4, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(h4, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(c32, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(h32, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(h3, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(c22, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(h22, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), c1, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), h1, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), c21, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), h21, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), c31, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), h31, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), c4, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), h4, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), c32, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), h32, cal(p['HC'][0], p['HA'][0], i), cal(p['HC'][1], p['HA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), h3, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), c22, cal(p['0_C'][0], p['CA'][0], i), cal(p['0_C'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), h22, cal(p['0_H'][0], p['HA'][0], i), cal(p['0_H'][1], p['HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', c1), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c1, h1), cal(p['CT_HC'][0], p['HC_sCA'][0], i), cal(p['CT_HC'][1], p['HC_sCA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c1, c21), cal(p['CA_sCT'][0], p['CA_CA'][0], i), cal(p['CA_sCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c1, c22), cal(p['CA_mCT'][0], p['CA_CA'][0], i), cal(p['CA_mCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c21, h21), cal(p['HA_sCT'][0], p['CA_HA'][0], i), cal(p['HA_sCT'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c21, c31), cal(p['CA_mCT'][0], p['CA_CA'][0], i), cal(p['CA_mCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c31, h31), cal(p['HA_mCT'][0], p['CA_HA'][0], i), cal(p['HA_mCT'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c31, c4), cal(p['CA_mCT'][0], p['CA_CA'][0], i), cal(p['CA_mCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c4, h4), cal(p['CT_HC'][0], p['HC_sCA'][0], i), cal(p['CT_HC'][1], p['HC_sCA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c4, 'CT'), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c4, c32), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c32, h32), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c32, h3), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c32, c22), cal(p['CA_mCT'][0], p['CA_CA'][0], i), cal(p['CA_mCT'][1], p['CA_CA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(c22, h22), cal(p['HA_mCT'][0], p['CA_HA'][0], i), cal(p['HA_mCT'][1], p['CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', c1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', c1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', 'CT', c1), cal(p['C_C_S'][0], p['CA_C_S'][0], i), cal(p['C_C_S'][1], p['CA_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', c1, h1), cal(p['C_C_H'][0], p['F_CT_CA_CA'][0], i), cal(p['C_C_H'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', c1, c21), cal(p['C_C_C'][0], p['F_CT_CA_CA'][0], i), cal(p['C_C_C'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', c1, c22), cal(p['C_C_C'][0], p['F_CT_CA_CA'][0], i), cal(p['C_C_C'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h1, c1, c21), cal(p['C_C_H'][0], p['Close'][0], i), cal(p['C_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h1, c1, h1), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h1, c1, c22), cal(p['C_C_H'][0], p['F_CA_CA_CA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c21, c1, c22), cal(p['Close'][0], p['F_CA_CA_CA'][0], i), cal(p['Close'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c1, c21, h21), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c1, c21, c31), cal(p['C_C_C'][0], p['F_CA_CA_CA'][0], i), cal(p['C_C_C'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h21, c21, c31), cal(p['C_C_C'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_C'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c21, c31, c4), cal(p['Dritt'][0], p['F_CA_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c21, c31, h31), cal(p['Dritt'][0], p['F_CA_CA_HA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h31, c31, c4), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c32, c4, h4), cal(p['C_C_H'][0], p['Close'][0], i), cal(p['C_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h4, c4, h4), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c31, c4, c32), cal(p['Close'][0], p['F_CA_CA_CA'][0], i), cal(p['Close'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h4, c4, c31), cal(p['C_C_H'][0], p['F_CA_CA_CA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c4, c32, h32), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c4, c32, h3), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c4, c32, c22), cal(p['C_C_C'][0], p['F_CA_CA_CA'][0], i), cal(p['C_C_C'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h3, c32, c22), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h32, c32, c22), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(h32, c32, h3), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c32, c22, c1), cal(p['Dritt'][0], p['F_CA_CA_CA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c32, c22, h22), cal(p['Close'][0], p['F_CA_CA_HA'][0], i), cal(p['Close'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(c1, c22, h22), cal(p['Dritt'][0], p['F_CA_CA_HA'][0], i), cal(p['Dritt'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', c4), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', c4), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', 'CT', c4), cal(p['C_C_S'][0], p['CA_C_S'][0], i), cal(p['C_C_S'][1], p['CA_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', c4, h4), cal(p['C_C_H'][0], p['F_CT_CA_CA'][0], i), cal(p['C_C_H'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', c4, c31), cal(p['C_C_C'][0], p['F_CT_CA_CA'][0], i), cal(p['C_C_C'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', c4, c32), cal(p['C_C_C'][0], p['F_CT_CA_CA'][0], i), cal(p['C_C_C'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', c1, h1), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', c1, c21), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', c1, c22), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', c1, h1), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', c1, c21), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', c1, c22), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', c1, h1), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', c1, c21), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', c1, c22), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', c1, c21, h21), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', c1, c22, h22), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', c1, c21, c31), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', c1, c22, c32), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h1, c1, c21, h21), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h1, c1, c22, h22), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h1, c1, c21, c31), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h1, c1, c22, c32), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c1, c21, c31, h31), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c1, c21, c31, c4), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h21, c21, c31, h31), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h21, c21, c31, c4), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c21, c31, c4, c32), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c21, c31, c4, 'CT'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h31, c31, c4, c32), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h31, c31, c4, 'CT'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h32, c32, c4, 'CT'), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h3, c32, c4, 'CT'), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c21, c31, c4, h4), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h31, c31, c4, h4), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c22, c32, c4, h4), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h32, c32, c4, h4), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h3, c32, c4, h4), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h3, c32, c4, c31), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h3, c32, c22, h22), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h3, c32, c22, c1), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', c4, c32, c22), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h32, c32, c22, h22), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(h32, c32, c22, c1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c4, c32, c22, c1), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c4, c32, c22, h22), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c22, c1, c21, h21), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c22, c1, c21, c31), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c21, c1, c22, h22), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c21, c1, c22, c32), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c31, c4, c32, h32), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(c31, c4, c32, c22), cal(p['Ring_0'][0], p['Ring_Dihe'][0], i), cal(p['Ring_0'][1], p['Ring_Dihe'][1], i), cal(p['Ring_0'][2], p['Ring_Dihe'][2], i), cal(p['Ring_0'][3], p['Ring_Dihe'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', c4, h4), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', c4, c31), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', c4, c32), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', c4, h4), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', c4, c31), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', c4, c32), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', c4, h4), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', c4, c31), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', 'CT', c4, c32), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), c1, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), h1, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), c21, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), h21, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), c31, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), h31, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), c4, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), h4, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), c32, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), h32, cal(p['HC'][2], p['HA'][2], i), cal(p['HC'][3], p['HA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), h3, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), c22, cal(p['0_C'][2], p['CA'][2], i), cal(p['0_C'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), h22, cal(p['0_H'][2], p['HA'][2], i), cal(p['0_H'][3], p['HA'][3], i))
