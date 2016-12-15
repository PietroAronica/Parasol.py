# VAL to GLN Mutation

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
        with open('Param_files/AminoAcid/GLN.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		i = float(i)
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HB2 :{}@HG21 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HG1 :{}@OE1 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HG1 :{}@NE2 0 0'.format(vxi, vxi)).execute()
		changeLJPair(parm, ':{}@HG1 :{}@HE21 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), fc['HB2']/10*i*i/10).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB']+((fc['HB3']-bc['HB'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG2'.format(vxi), bc['CG2']-(bc['CG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HG21'.format(vxi), bc['HG21']-(bc['HG21']/10)*i).execute()
                change(parm, 'charge', ':{}@HG22'.format(vxi), bc['HG22']-(bc['HG22']/10)*i).execute()
                change(parm, 'charge', ':{}@HG23'.format(vxi), bc['HG23']-(bc['HG23']/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG1']+((fc['CG']-bc['CG1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vxi), bc['HG11']-(bc['HG11']/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), bc['HG12']+((fc['HG2']-bc['HG12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), bc['HG13']+((fc['HG3']-bc['HG13'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), fc['CD']/10*i*i/10).execute()
                change(parm, 'charge', ':{}@OE1'.format(vxi), fc['OE1']/10*i*i/10).execute()
                change(parm, 'charge', ':{}@NE2'.format(vxi), fc['NE2']/10*i*i/10).execute()
                change(parm, 'charge', ':{}@HE21'.format(vxi), fc['HE21']/10*i*i/10).execute()
                change(parm, 'charge', ':{}@HE22'.format(vxi), fc['HE22']/10*i*i/10).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		#print printDetails(parm, ':VXI')
		#d = netCharge(parm).execute()
		#change(parm, 'charge', ':PN', '{:.3f}'.format(-d)).execute() 
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CG1 = struct.residue_dict[aa].atom_dict['CG1']
        HG11 = struct.residue_dict[aa].atom_dict['HG11']
        CG2 = struct.residue_dict[aa].atom_dict['CG2']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'CB' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.superimposed1('HB2', CG2))
			elif atom.get_name() == 'HB' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HB3'))
			elif atom.get_name() == 'CG1' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('CG'))
			elif atom.get_name() == 'HG11' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HG1'))
			elif atom.get_name() == 'HG12' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HG2'))
			elif atom.get_name() == 'HG13' and res.get_resname() == vxi:
                        	pdb.write(atom.change_name('HG3'))
                        	pdb.write(atom.halfway_between('CD', CG1, HG11))
                        	pdb.write(atom.superimposed1('OE1', HG11))
                        	pdb.write(atom.superimposed2('NE2', HG11))
                                pdb.write(atom.halfway_after1('HE21', CG1, HG11))
                                pdb.write(atom.halfway_after2('HE22', CG1, HG11))	
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

def lib_make(ff, outputfile, vxi='VXI', metcar='dc', methyd='dh', hydhyd1='mh', amicar='ac', amioxy='ao', aminit='an', amihyd='ah', hydhyd2='sh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/VAL-GLN.pdb\n"%vxi)
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
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "O"\n'%vxi)
	ctrl.write('set %s.1.18 element "N"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "C"\n'%vxi)
	ctrl.write('set %s.1.22 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG2"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG21"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG22"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG23"\n'%vxi)
	ctrl.write('set %s.1.12 name "CG"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.14 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.15 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.16 name "CD"\n'%vxi)
	ctrl.write('set %s.1.17 name "OE1"\n'%vxi)
	ctrl.write('set %s.1.18 name "NE2"\n'%vxi)
	ctrl.write('set %s.1.19 name "HE21"\n'%vxi)
	ctrl.write('set %s.1.20 name "HE22"\n'%vxi)
	ctrl.write('set %s.1.21 name "C"\n'%vxi)
	ctrl.write('set %s.1.22 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, hydhyd1))
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.12 type "CT"\n'%vxi)
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, hydhyd2))
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "HC"\n'%vxi)
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, amicar))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, amioxy))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, aminit))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, amihyd))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, amihyd))
	ctrl.write('set %s.1.21 type "C"\n'%vxi)
	ctrl.write('set %s.1.22 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.22\n'%(vxi, vxi))
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

def stock_add_to_all(metcar='dc', methyd='dh', hydhyd1='mh', amicar='ac', amioxy='ao', aminit='an', amihyd='ah', hydhyd2='sh'):
	Frcmod_creator.make_hyb()
        Frcmod_creator.TYPE_insert(metcar, 'C', 'sp3')
        Frcmod_creator.TYPE_insert(methyd, 'H', 'sp3')
        Frcmod_creator.TYPE_insert(hydhyd1, 'H', 'sp3')	
        Frcmod_creator.TYPE_insert(amicar, 'C', 'sp2')
        Frcmod_creator.TYPE_insert(amioxy, 'O', 'sp2')
        Frcmod_creator.TYPE_insert(aminit, 'N', 'sp2')
        Frcmod_creator.TYPE_insert(amihyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd2, 'H', 'sp3')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), cal(p['CT_CT'][0], p['CT_mH'][0], i), cal(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd1), cal(p['HC_sC'][0], p['CT_HC'][0], i), cal(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), cal(p['CT_HC'][0], p['HC_mH'][0], i), cal(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), cal(p['C_C_H'][0], p['Dritt'][0], i), cal(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', metcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', hydhyd1), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', metcar), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', metcar, methyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', metcar, methyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, cal(p['CT'][2], p['0_C'][2], i), cal(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))

                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), amicar, lac(p['C'][0], p['0_C'][0], i), lac(p['C'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), amioxy, lac(p['O2'][0], p['0_O'][0], i), lac(p['O2'][1], p['0_O'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), aminit, lac(p['NA'][0], p['0_N'][0], i), lac(p['NA'][1], p['0_N'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), amihyd, lac(p['H'][0], p['0_H'][0], i), lac(p['H'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, lac(p['0_H'][0], p['HC'][0], i), lac(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', amicar), lac(p['CT_C'][0], p['CT_mH'][0], i), lac(p['CT_C'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd2), lac(p['HC_sC2'][0], p['CT_HC'][0], i), lac(p['HC_sC2'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(amicar, amioxy), lac(p['C_O'][0], p['O_mH'][0], i), lac(p['C_O'][1], p['O_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(amicar, aminit), lac(p['C_N'][0], p['N_mH'][0], i), lac(p['C_N'][1], p['N_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(aminit, amihyd), lac(p['NA_H'][0], p['H_mHC'][0], i), lac(p['NA_H'][1], p['H_mHC'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', amicar, amioxy), lac(p['C_C_O'][0], p['Dritt'][0], i), lac(p['C_C_O'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', amicar, aminit), lac(p['C_C_N'][0], p['Dritt'][0], i), lac(p['C_C_N'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amioxy, amicar, aminit), lac(p['O_C_N'][0], p['Close'][0], i), lac(p['O_C_N'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', amicar), lac(p['CT_CT_C'][0], p['C_C_H'][0], i), lac(p['CT_CT_C'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', amicar), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', hydhyd2), lac(p['H_C_H'][0], p['H_C_H'][0], i), lac(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', amicar), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd2), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amicar, aminit, amihyd), lac(p['F_CA_CA_HA'][0], p['Dritt'][0], i), lac(p['F_CA_CA_HA'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amihyd, aminit, amihyd), lac(p['H_N_H'][0], p['Close'][0], i), lac(p['H_N_H'][1], p['Close'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, amioxy), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), lac(p['C_C_C_N_1'][0], p['0_3'][0], i), lac(p['C_C_C_N_1'][1], p['0_3'][1], i), lac(p['C_C_C_N_1'][2], p['0_3'][2], i), lac(p['C_C_C_N_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), lac(p['C_C_C_N_2'][0], p['0_8'][0], i), lac(p['C_C_C_N_2'][1], p['0_8'][1], i), lac(p['C_C_C_N_2'][2], p['0_8'][2], i), lac(p['C_C_C_N_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), lac(p['C_C_C_N_3'][0], p['0_2'][0], i), lac(p['C_C_C_N_3'][1], p['0_2'][1], i), lac(p['C_C_C_N_3'][2], p['0_2'][2], i), lac(p['C_C_C_N_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), lac(p['C_C_C_N_4'][0], p['0_7'][0], i), lac(p['C_C_C_N_4'][1], p['0_7'][1], i), lac(p['C_C_C_N_4'][2], p['0_7'][2], i), lac(p['C_C_C_N_4'][3], p['0_7'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', amicar, aminit), lac(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), lac(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), lac(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), lac(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', amicar, amioxy), lac(p['H_C_C_O_1'][0], p['0_10'][0], i), lac(p['H_C_C_O_1'][1], p['0_10'][1], i), lac(p['H_C_C_O_1'][2], p['0_10'][2], i), lac(p['H_C_C_O_1'][3], p['0_10'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', amicar, amioxy), lac(p['H_C_C_O_2'][0], p['0_8'][0], i), lac(p['H_C_C_O_2'][1], p['0_8'][1], i), lac(p['H_C_C_O_2'][2], p['0_8'][2], i), lac(p['H_C_C_O_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', amicar, amioxy), lac(p['H_C_C_O_3'][0], p['0_9'][0], i), lac(p['H_C_C_O_3'][1], p['0_9'][1], i), lac(p['H_C_C_O_3'][2], p['0_9'][2], i), lac(p['H_C_C_O_3'][3], p['0_9'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', amicar, amioxy), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', amicar, aminit), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amihyd, aminit, amicar, amioxy), lac(p['O_C_N_H_1'][0], p['0_3'][0], i), lac(p['O_C_N_H_1'][1], p['0_3'][1], i), lac(p['O_C_N_H_1'][2], p['0_3'][2], i), lac(p['O_C_N_H_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amihyd, aminit, amicar, amioxy), lac(p['O_C_N_H_2'][0], p['0_11'][0], i), lac(p['O_C_N_H_2'][1], p['0_11'][1], i), lac(p['O_C_N_H_2'][2], p['0_11'][2], i), lac(p['O_C_N_H_2'][3], p['0_11'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amihyd, aminit, amicar, 'CT'), lac(p['X_C_N_X'][0], p['Ring_0'][0], i), lac(p['X_C_N_X'][1], p['Ring_0'][1], i), lac(p['X_C_N_X'][2], p['Ring_0'][2], i), lac(p['X_C_N_X'][3], p['Ring_0'][3], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', amicar, amioxy), lac(p['Car_imp'][0], p['Imp_0'][0], i), lac(p['Car_imp'][1], p['Imp_0'][1], i), lac(p['Car_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', aminit, amihyd), lac(p['Ami_imp'][0], p['Imp_0'][0], i), lac(p['Ami_imp'][1], p['Imp_0'][1], i), lac(p['Ami_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amicar, 'CT', aminit, amihyd), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amicar, 'CT', aminit, amioxy), lac(p['Ring_imp'][0], p['Imp_0'][0], i), lac(p['Ring_imp'][1], p['Imp_0'][1], i), lac(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), amicar, lac(p['C'][2], p['0_C'][2], i), lac(p['C'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), amioxy, lac(p['O2'][2], p['0_O'][2], i), lac(p['O2'][3], p['0_O'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), aminit, lac(p['NA'][2], p['0_N'][2], i), lac(p['NA'][3], p['0_N'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), amihyd, lac(p['H'][2], p['0_H'][2], i), lac(p['H'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, lac(p['0_H'][2], p['HC'][2], i), lac(p['0_H'][3], p['HC'][3], i))
