# X5X to XnX Mutation

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
        with open('Param_files/AminoAcid/XnX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HD11'.format(vxi), ':{}@CB1'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@HD21'.format(vxi), ':{}@CB2'.format(vxi), '0', '0').execute()
		changeLJPair(parm, ':{}@OE2'.format(vxi), ':{}@OE1'.format(vxi), '0', '0').execute()
                change(parm, 'charge', ':{}@CA1'.format(vxi), bc['C1']+((fc['CA1']-bc['C1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA12'.format(vxi), bc['H12']+((fc['HA12']-bc['H12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA13'.format(vxi), bc['H13']+((fc['HA13']-bc['H13'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB1'.format(vxi), bc['C2']+((fc['CB1']-bc['C2'])/10)*i).execute()
                change(parm, 'charge', ':{}@CC1'.format(vxi), bc['H22']+((fc['CC1']-bc['H22'])/10)*i).execute()
                change(parm, 'charge', ':{}@H23'.format(vxi), bc['H23']-(bc['H23']/10)*i).execute()
                change(parm, 'charge', ':{}@CD1'.format(vxi), (fc['CD1']/10)*i).execute()
                change(parm, 'charge', ':{}@HD11'.format(vxi), (fc['HD11']/10)*i).execute()
                change(parm, 'charge', ':{}@HD12'.format(vxi), (fc['HD12']/10)*i).execute()
                change(parm, 'charge', ':{}@HD13'.format(vxi), (fc['HD13']/10)*i).execute()
                change(parm, 'charge', ':{}@CE1'.format(vxi), (fc['CE1']/10)*i).execute()
                change(parm, 'charge', ':{}@OE1'.format(vxi), (fc['OE1']/10)*i).execute()
                change(parm, 'charge', ':{}@NF2'.format(vxi), bc['H32']+((fc['NF2']-bc['H32'])/10)*i).execute()
                change(parm, 'charge', ':{}@NF1'.format(vxi), bc['C3']+((fc['NF1']-bc['C3'])/10)*i).execute()
                change(parm, 'charge', ':{}@H33'.format(vxi), bc['H33']-(bc['H33']/10)*i).execute()
                change(parm, 'charge', ':{}@CE2'.format(vxi), (fc['CE2']/10)*i).execute()
                change(parm, 'charge', ':{}@OE2'.format(vxi), (fc['OE2']/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), (fc['CD2']/10)*i).execute()
                change(parm, 'charge', ':{}@HD21'.format(vxi), (fc['HD21']/10)*i).execute()
                change(parm, 'charge', ':{}@HD22'.format(vxi), (fc['HD22']/10)*i).execute()
                change(parm, 'charge', ':{}@HD23'.format(vxi), (fc['HD23']/10)*i).execute()
                change(parm, 'charge', ':{}@CC2'.format(vxi), bc['H42']+((fc['CC2']-bc['H42'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB2'.format(vxi), bc['C4']+((fc['CB2']-bc['C4'])/10)*i).execute()
                change(parm, 'charge', ':{}@H43'.format(vxi), bc['H43']-(bc['H43']/10)*i).execute()
                change(parm, 'charge', ':{}@CA2'.format(vxi), bc['C5']+((fc['CA2']-bc['C5'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA22'.format(vxi), bc['H52']+((fc['HA22']-bc['H52'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA23'.format(vxi), bc['H53']+((fc['HA23']-bc['H53'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        C2 = struct.residue_dict[aa].atom_dict['C2']
        H22 = struct.residue_dict[aa].atom_dict['H22']
        H32 = struct.residue_dict[aa].atom_dict['H32']
        C4 = struct.residue_dict[aa].atom_dict['C4']
        H42 = struct.residue_dict[aa].atom_dict['H42']
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
                        	pdb.write(atom.change_name('CB1'))
			elif atom.get_name() == 'H22' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CC1'))
			elif atom.get_name() == 'H23' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('CD1', C2, H22))
                        	pdb.write(atom.superimposed1('HD11', C2))
                        	pdb.write(atom.superimposed2('HD12', C2))
                        	pdb.write(atom.superimposed3('HD13', C2))
			elif atom.get_name() == 'C3' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('NF1'))
                        	pdb.write(atom.halfway_between('CE1', H22, H32))
                        	pdb.write(atom.superimposed1('OE1', H32))
			elif atom.get_name() == 'H32' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('NF2'))
			elif atom.get_name() == 'H33' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('CE2', H42, H32))
                        	pdb.write(atom.superimposed2('OE2', H32))
			elif atom.get_name() == 'H42' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CC2'))
			elif atom.get_name() == 'H43' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.halfway_between('CD2', C4, H42))
                        	pdb.write(atom.superimposed1('HD21', C4))
                        	pdb.write(atom.superimposed2('HD22', C4))
                        	pdb.write(atom.superimposed3('HD23', C4))
			elif atom.get_name() == 'C4' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CB2'))
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
 	cb = var[0]
 	hb = var[1]
 	cc = var[2]
        cd = var[3]
        hd = var[4]
        ce = var[5]
        oe = var[6]
        n1 = var[7]
        n2 = var[8]
        hn = var[9]
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/X5X-XnX.pdb\n"%vxi)
	ctrl.write('set %s.1.1  element "C"\n'%vxi)
	ctrl.write('set %s.1.2  element "H"\n'%vxi)
	ctrl.write('set %s.1.3  element "H"\n'%vxi)
	ctrl.write('set %s.1.4  element "C"\n'%vxi)
	ctrl.write('set %s.1.5  element "H"\n'%vxi)
	ctrl.write('set %s.1.6  element "C"\n'%vxi)
	ctrl.write('set %s.1.7  element "C"\n'%vxi)
	ctrl.write('set %s.1.8  element "H"\n'%vxi)
	ctrl.write('set %s.1.9  element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "O"\n'%vxi)
	ctrl.write('set %s.1.13 element "N"\n'%vxi)
	ctrl.write('set %s.1.14 element "N"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "O"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "C"\n'%vxi)
	ctrl.write('set %s.1.24 element "H"\n'%vxi)
	ctrl.write('set %s.1.25 element "C"\n'%vxi)
	ctrl.write('set %s.1.26 element "H"\n'%vxi)
	ctrl.write('set %s.1.27 element "H"\n'%vxi)
	ctrl.write('set %s.1.1  name "CA1"\n'%vxi)
	ctrl.write('set %s.1.2  name "HA12"\n'%vxi)
	ctrl.write('set %s.1.3  name "HA13"\n'%vxi)
	ctrl.write('set %s.1.4  name "CB1"\n'%vxi)
	ctrl.write('set %s.1.5  name "H23"\n'%vxi)
	ctrl.write('set %s.1.6  name "CC1"\n'%vxi)
	ctrl.write('set %s.1.7  name "CD1"\n'%vxi)
	ctrl.write('set %s.1.8  name "HD11"\n'%vxi)
	ctrl.write('set %s.1.9  name "HD12"\n'%vxi)
	ctrl.write('set %s.1.10 name "HD13"\n'%vxi)
	ctrl.write('set %s.1.11 name "CE1"\n'%vxi)
	ctrl.write('set %s.1.12 name "OE1"\n'%vxi)
	ctrl.write('set %s.1.13 name "NF1"\n'%vxi)
	ctrl.write('set %s.1.14 name "NF2"\n'%vxi)
	ctrl.write('set %s.1.15 name "H33"\n'%vxi)
	ctrl.write('set %s.1.16 name "CE2"\n'%vxi)
	ctrl.write('set %s.1.17 name "OE2"\n'%vxi)
	ctrl.write('set %s.1.18 name "CC2"\n'%vxi)
	ctrl.write('set %s.1.19 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.20 name "HD21"\n'%vxi)
	ctrl.write('set %s.1.21 name "HD22"\n'%vxi)
	ctrl.write('set %s.1.22 name "HD23"\n'%vxi)
	ctrl.write('set %s.1.23 name "CB2"\n'%vxi)
	ctrl.write('set %s.1.24 name "H43"\n'%vxi)
	ctrl.write('set %s.1.25 name "CA2"\n'%vxi)
	ctrl.write('set %s.1.26 name "HA22"\n'%vxi)
	ctrl.write('set %s.1.27 name "HA23"\n'%vxi)
	ctrl.write('set %s.1.1 type "CT"\n'%vxi)
	ctrl.write('set %s.1.2 type "HC"\n'%vxi)
	ctrl.write('set %s.1.3 type "HC"\n'%vxi)
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, cb))
	ctrl.write('set %s.1.5 type "%s"\n'%(vxi, hb))
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, cc))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, cd))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, hd))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, hd))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, hd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, ce))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, oe))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, n1))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, n2))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, hn))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, ce))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, oe))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, cc))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, cd))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, hd))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, hd))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, hd))
	ctrl.write('set %s.1.23 type "%s"\n'%(vxi, cb))
	ctrl.write('set %s.1.24 type "%s"\n'%(vxi, hb))
	ctrl.write('set %s.1.25 type "CT"\n'%vxi)
	ctrl.write('set %s.1.26 type "HC"\n'%vxi)
	ctrl.write('set %s.1.27 type "HC"\n'%vxi)
	ctrl.write('bond %s.1.1  %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1  %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1  %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4  %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4  %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.4  %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6  %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.6  %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7  %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7  %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7  %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.23 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.26\n'%(vxi, vxi))
	ctrl.write('bond %s.1.25 %s.1.27\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.CA1\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.CA2\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.CA1\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.CA2\n'%(vxi, vxi))
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
 	cb = var[0]
 	hb = var[1]
 	cc = var[2]
        cd = var[3]
        hd = var[4]
        ce = var[5]
        oe = var[6]
        n1 = var[7]
        n2 = var[8]
        hn = var[9]
        Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(cb, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(hb, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(cc, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(cd, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(hd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(ce, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(oe, 'O', 'sp2')
        Frcmod_creator.TYPE_insert(n1, 'N', 'sp2')
        Frcmod_creator.TYPE_insert(n2, 'N', 'sp2')
        Frcmod_creator.TYPE_insert(hn, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cb, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cc, cal(p['HC'][0], p['CA'][0], i), cal(p['HC'][1], p['CA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hb, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), cd, cal(p['0_C'][0], p['CT'][0], i), cal(p['0_C'][1], p['CT'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hd, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), ce, cal(p['0_C'][0], p['C'][0], i), cal(p['0_C'][1], p['C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), oe, cal(p['0_O'][0], p['O2'][0], i), cal(p['0_O'][1], p['O2'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), n1, cal(p['CT'][0], p['NA'][0], i), cal(p['CT'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), n2, cal(p['HC'][0], p['NA'][0], i), cal(p['HC'][1], p['NA'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hn, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', cb), cal(p['CT_CT'][0], p['CT_CC'][0], i), cal(p['CT_CT'][1], p['CT_CC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cb, cc), cal(p['CT_HC'][0], p['CC_CD'][0], i), cal(p['CT_HC'][1], p['CC_CD'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cb, n1), cal(p['CT_CT'][0], p['CC_NA'][0], i), cal(p['CT_CT'][1], p['CC_NA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cb, hb), cal(p['CT_HC'][0], p['CC_sHC'][0], i), cal(p['CT_HC'][1], p['CC_sHC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cc, cd), cal(p['CC_mHC'][0], p['CT_CC'][0], i), cal(p['CC_mHC'][1], p['CT_CC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cd, hd), cal(p['CC_mHC'][0], p['CT_HC'][0], i), cal(p['CC_mHC'][1], p['CT_HC'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(cc, ce), cal(p['0_0'][0], p['CC_C'][0], i), cal(p['0_0'][1], p['CC_C'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce, oe), cal(p['0_0'][0], p['C_O'][0], i), cal(p['0_0'][1], p['C_O'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(ce, n2), cal(p['0_0'][0], p['C_NA'][0], i), cal(p['0_0'][1], p['C_NA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(n1, n2), cal(p['CT_HC'][0], p['N_NA'][0], i), cal(p['CT_HC'][1], p['N_NA'][1], i))
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(n1, hn), cal(p['CT_HC'][0], p['N_sHC'][0], i), cal(p['CT_HC'][1], p['N_sHC'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', cb), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', cb), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', 'CT', cb), cal(p['C_C_S'][0], p['S_CT_CC'][0], i), cal(p['C_C_S'][1], p['S_CT_CC'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', cb, n1), cal(p['C_C_C'][0], p['F_CT_CA_CA'][0], i), cal(p['C_C_C'][1], p['F_CT_CA_CA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', cb, hb), cal(p['C_C_H'][0], p['CT_CC_CD'][0], i), cal(p['C_C_H'][1], p['CT_CC_CD'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', cb, cc), cal(p['C_C_H'][0], p['CT_CC_CD'][0], i), cal(p['C_C_H'][1], p['CT_CC_CD'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hb, cb, cc), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(n1, cb, cc), cal(p['C_C_H'][0], p['CD_CC_NA'][0], i), cal(p['C_C_H'][1], p['CD_CC_NA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(n1, cb, hb), cal(p['C_C_H'][0], p['CD_CC_NA'][0], i), cal(p['C_C_H'][1], p['CD_CC_NA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hn, n1, cb), cal(p['C_C_H'][0], p['CC_NA_N'][0], i), cal(p['C_C_H'][1], p['CC_NA_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(n2, n1, cb), cal(p['C_C_H'][0], p['CC_NA_N'][0], i), cal(p['C_C_H'][1], p['CC_NA_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cb, n1, cb), cal(p['C_C_C'][0], p['CC_NA_CC'][0], i), cal(p['C_C_C'][1], p['CC_NA_CC'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cb, cc, cd), cal(p['Close'][0], p['CT_CC_CD'][0], i), cal(p['Close'][1], p['CT_CC_CD'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cc, cd, hd), cal(p['Dritt'][0], p['C_C_H'][0], i), cal(p['Dritt'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hd, cd, hd), cal(p['Close'][0], p['H_C_H'][0], i), cal(p['Dritt'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cb, cc, ce), cal(p['180_Close'][0], p['C_CD_CC'][0], i), cal(p['180_Close'][1], p['C_CD_CC'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cd, cc, ce), cal(p['180_Close'][0], p['CT_CD_C'][0], i), cal(p['180_Close'][1], p['CT_CD_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cc, ce, oe), cal(p['180_Close'][0], p['CD_C_O'][0], i), cal(p['180_Close'][1], p['CD_C_O'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(n2, ce, oe), cal(p['180_Close'][0], p['O_C_N'][0], i), cal(p['180_Close'][1], p['O_C_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(cc, ce, n2), cal(p['180_Close'][0], p['CD_C_N'][0], i), cal(p['180_Close'][1], p['CD_C_N'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce, n2, n1), cal(p['180_Close'][0], p['C_N_NA'][0], i), cal(p['180_Close'][1], p['C_N_NA'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(ce, n2, ce), cal(p['180_Close'][0], p['C_N_C'][0], i), cal(p['180_Close'][1], p['C_N_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(n2, n1, hn), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'CT', cb, 'X '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', cc, cd, 'X '), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', cb, cc, 'X '), cal(p['Ring_0'][0], p['X_CB_CC_X'][0], i), cal(p['Ring_0'][1], p['X_CB_CC_X'][1], i), cal(p['Ring_0'][2], p['X_CB_CC_X'][2], i), cal(p['Ring_0'][3], p['X_CB_CC_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', cb, n1, 'X '), cal(p['Ring_0'][0], p['X_CC_NA_X'][0], i), cal(p['Ring_0'][1], p['X_CC_NA_X'][1], i), cal(p['Ring_0'][2], p['X_CC_NA_X'][2], i), cal(p['Ring_0'][3], p['X_CC_NA_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', cc, ce, 'X '), cal(p['Ring_0'][0], p['X_CC_C_X'][0], i), cal(p['Ring_0'][1], p['X_CC_C_X'][1], i), cal(p['Ring_0'][2], p['X_CC_C_X'][2], i), cal(p['Ring_0'][3], p['X_CC_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', n1, n2, 'X '), cal(p['0_19'][0], p['X_NA_N_X'][0], i), cal(p['0_19'][1], p['X_NA_N_X'][1], i), cal(p['0_19'][2], p['X_NA_N_X'][2], i), cal(p['0_19'][3], p['X_NA_N_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', ce, n2, 'X '), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cb, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cc, cal(p['HC'][2], p['CA'][2], i), cal(p['HC'][3], p['CA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hb, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), cd, cal(p['0_C'][2], p['CT'][2], i), cal(p['0_C'][3], p['CT'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hd, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), ce, cal(p['0_C'][2], p['C'][2], i), cal(p['0_C'][3], p['C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), oe, cal(p['0_O'][2], p['O2'][2], i), cal(p['0_O'][3], p['O2'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), n1, cal(p['CT'][2], p['NA'][2], i), cal(p['CT'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), n2, cal(p['HC'][2], p['NA'][2], i), cal(p['HC'][3], p['NA'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hn, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
