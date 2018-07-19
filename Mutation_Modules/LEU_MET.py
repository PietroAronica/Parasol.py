# LEU to MET Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command(vxi='VXI', lipid='No'):
	bc = {}
        with open('Param_files/AminoAcid/LEU.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/MET.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                changeLJPair(parm, ':{}@HG3 :{}@HD21 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB2']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG']+((fc['HG2']-bc['HG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), fc['HG3']/10*i).execute()
                change(parm, 'charge', ':{}@SD'.format(vxi), fc['SD']/10*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']-(bc['CD2']/10)*i).execute()
                change(parm, 'charge', ':{}@HD21'.format(vxi), bc['HD21']-(bc['HD21']/10)*i).execute()
                change(parm, 'charge', ':{}@HD22'.format(vxi), bc['HD22']-(bc['HD22']/10)*i).execute()
                change(parm, 'charge', ':{}@HD23'.format(vxi), bc['HD23']-(bc['HD23']/10)*i).execute()
                change(parm, 'charge', ':{}@CE'.format(vxi), bc['CD1']+((fc['CE']-bc['CD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE1'.format(vxi), bc['HD11']+((fc['HE1']-bc['HD11'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vxi), bc['HD12']+((fc['HE2']-bc['HD12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vxi), bc['HD13']+((fc['HE3']-bc['HD13'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CG = struct.residue_dict[aa].atom_dict['CG']
        CD1 = struct.residue_dict[aa].atom_dict['CD1']
        CD2 = struct.residue_dict[aa].atom_dict['CD2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HG' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HG2'))
                        	pdb.write(atom.superimposed1('HG3', CD2))
                        	pdb.write(atom.halfway_between('SD', CG, CD1))
			elif atom.get_name() == 'CD1' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CE'))
			elif atom.get_name() == 'HD11' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HE1'))
			elif atom.get_name() == 'HD12' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HE2'))
			elif atom.get_name() == 'HD13' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HE3'))
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
	metcar = var[0]
	methyd = var[1]
	hydhyd = var[2]
	estsul = var[3]
	mthhyd = var[4]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/LEU-MET.pdb\n"%vxi)
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
	ctrl.write('set %s.1.11 element "S"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "C"\n'%vxi)
	ctrl.write('set %s.1.21 element "O"\n'%vxi)
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
	ctrl.write('set %s.1.11 name "SD"\n'%vxi)
	ctrl.write('set %s.1.12 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.13 name "HD21"\n'%vxi)
	ctrl.write('set %s.1.14 name "HD22"\n'%vxi)
	ctrl.write('set %s.1.15 name "HD23"\n'%vxi)
	ctrl.write('set %s.1.16 name "CE"\n'%vxi)
	ctrl.write('set %s.1.17 name "HE1"\n'%vxi)
	ctrl.write('set %s.1.18 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.19 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.20 name "C"\n'%vxi)
	ctrl.write('set %s.1.21 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "CT"\n'%vxi)
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, mthhyd))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, hydhyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, estsul))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.16 type "CT"\n'%vxi)
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, mthhyd))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, mthhyd))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, mthhyd))
	ctrl.write('set %s.1.20 type "C"\n'%vxi)
	ctrl.write('set %s.1.21 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.20 %s.1.21\n'%(vxi, vxi))
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
	metcar = var[0]
	methyd = var[1]
	hydhyd = var[2]
	estsul = var[3]
	mthhyd = var[4]
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(estsul, 'S', 'sp3')
	Frcmod_creator.TYPE_insert(mthhyd, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), estsul, cal(p['0_S'][0], p['SH'][0], i), cal(p['0_S'][1], p['SH'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), mthhyd, cal(p['HC'][0], p['H1'][0], i), cal(p['HC'][1], p['H1'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', estsul), cal(p['SH_mCT'][0], p['CT_SH'][0], i), cal(p['SH_mCT'][1], p['CT_SH'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', mthhyd), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', estsul, 'CT'), cal(p['Dritt'][0], p['C_S_C'][0], i), cal(p['Dritt'][1], p['C_S_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(mthhyd, 'CT', estsul), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(estsul, 'CT', 'CT'), cal(p['C_C_C'][0], p['C_C_S'][0], i), cal(p['C_C_C'][1], p['C_C_S'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(mthhyd, 'CT', 'CT'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(mthhyd, 'CT', mthhyd), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', estsul, 'CT', mthhyd), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', estsul, 'CT', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), estsul, cal(p['0_S'][2], p['SH'][2], i), cal(p['0_S'][3], p['SH'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), mthhyd, cal(p['HC'][2], p['H1'][2], i), cal(p['HC'][3], p['H1'][3], i))

		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', estsul), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(metcar, 'CT', estsul), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(mthhyd, 'CT', hydhyd), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(mthhyd, 'CT', metcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', estsul, 'CT', hydhyd), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(mthhyd, 'CT', metcar, methyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(estsul, 'CT', metcar, methyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))

                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, cal(p['CT'][0], p['0_C'][0], i), cal(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), cal(p['CT_CT'][0], p['CT_mH'][0], i), cal(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), cal(p['HC_sC'][0], p['CT_HC'][0], i), cal(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), cal(p['CT_HC'][0], p['HC_mH'][0], i), cal(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), cal(p['C_C_H'][0], p['Dritt'][0], i), cal(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', metcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', hydhyd), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', metcar), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', metcar, methyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', estsul, 'CT', metcar), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, cal(p['CT'][2], p['0_C'][2], i), cal(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, cal(p['0_H'][2], p['H1'][2], i), cal(p['0_H'][3], p['H1'][3], i))
