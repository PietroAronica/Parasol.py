# ASN to THR Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
	bc = {}
        with open('Param_files/AminoAcid/ASN.param', 'r') as b:
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
                changeLJPair(parm, ':{}@HB :{}@OD1 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HB :{}@ND2 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HB :{}@HD21 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HB2 :{}@HG21 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HB3 :{}@HG1 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB'.format(vxi), (fc['HB']/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']-(bc['HB2']/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']-(bc['HB3']/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), (bc['CG']-(bc['CG']/10)*i)*(10-i)/10).execute()
                change(parm, 'charge', ':{}@OD1'.format(vxi), (bc['OD1']-(bc['OD1']/10)*i)*(10-i)/10).execute()
                change(parm, 'charge', ':{}@ND2'.format(vxi), (bc['ND2']-(bc['ND2']/10)*i)*(10-i)/10).execute()
                change(parm, 'charge', ':{}@HD21'.format(vxi), (bc['HD21']-(bc['HD21']/10)*i)*(10-i)/10).execute()
                change(parm, 'charge', ':{}@HD22'.format(vxi), (bc['HD22']-(bc['HD22']/10)*i)*(10-i)/10).execute()
                change(parm, 'charge', ':{}@CG2'.format(vxi), (fc['CG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HG21'.format(vxi), (fc['HG21']/10)*i).execute()
                change(parm, 'charge', ':{}@HG22'.format(vxi), (fc['HG22']/10)*i).execute()
                change(parm, 'charge', ':{}@HG23'.format(vxi), (fc['HG23']/10)*i).execute()
                change(parm, 'charge', ':{}@OG1'.format(vxi), (fc['OG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vxi), (fc['HG1']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CB = struct.residue_dict[aa].atom_dict['CB']
        HB2 = struct.residue_dict[aa].atom_dict['HB2']
        HB3 = struct.residue_dict[aa].atom_dict['HB3']
        CG = struct.residue_dict[aa].atom_dict['CG']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'CB' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.superimposed1('HB', CG))
			elif atom.get_name() == 'HD22' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                                pdb.write(atom.halfway_between('CG2', CB, HB2))
                                pdb.write(atom.superimposed1('HG21', HB2))
                                pdb.write(atom.superimposed2('HG22', HB2))
                                pdb.write(atom.superimposed3('HG23', HB2))
                                pdb.write(atom.halfway_between('OG1', CB, HB3))
                                pdb.write(atom.superimposed1('HG1', HB3))
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

def lib_make(ff, outputfile, vxi='VXI', metcar='dc', methyd='dh', hydhyd1='mh', alcoxy='ho', alchyd='hh', hydhyd2='sh', amicar='ac', amioxy='ao', aminit='an', amihyd='ah', hydhyd3='yh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/THR-ASN.pdb\n"%vxi)
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
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "O"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "C"\n'%vxi)
	ctrl.write('set %s.1.21 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.9 name "CG"\n'%vxi)
	ctrl.write('set %s.1.10 name "OD1"\n'%vxi)
	ctrl.write('set %s.1.11 name "ND2"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD21"\n'%vxi)
	ctrl.write('set %s.1.13 name "HD22"\n'%vxi)
	ctrl.write('set %s.1.14 name "CG2"\n'%vxi)
	ctrl.write('set %s.1.15 name "HG21"\n'%vxi)
	ctrl.write('set %s.1.16 name "HG22"\n'%vxi)
	ctrl.write('set %s.1.17 name "HG23"\n'%vxi)
	ctrl.write('set %s.1.18 name "OG1"\n'%vxi)
	ctrl.write('set %s.1.19 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.20 name "C"\n'%vxi)
	ctrl.write('set %s.1.21 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, hydhyd3))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, hydhyd1))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, hydhyd2))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, amicar))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, amioxy))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, aminit))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, amihyd))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, amihyd))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, alcoxy))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, alchyd))
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
	ctrl.write('bond %s.1.5 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
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

def lac(y, x, i):
        num = x+((y-x)/10)*i
        return num

def stock_add_to_all(metcar='dc', methyd='dh', hydhyd1='mh', alcoxy='ho', alchyd='hh', hydhyd2='sh', amicar='ac', amioxy='ao', aminit='an', amihyd='ah', hydhyd3='yh'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(metcar, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(methyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(alcoxy, 'O', 'sp3')
	Frcmod_creator.TYPE_insert(alchyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(amicar, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(amioxy, 'O', 'sp2')
	Frcmod_creator.TYPE_insert(aminit, 'N', 'sp2')
	Frcmod_creator.TYPE_insert(amihyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd3, 'H', 'sp3')
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
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][0], p['0_C'][0], i), lac(p['CT'][1], p['0_C'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, lac(p['0_H'][0], p['HC'][0], i), lac(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), lac(p['CT_CT'][0], p['CT_mH'][0], i), lac(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd1), lac(p['HC_sC'][0], p['CT_HC'][0], i), lac(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), lac(p['CT_HC'][0], p['HC_mH'][0], i), lac(p['CT_HC'][1], p['HC_mH'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', metcar), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), lac(p['C_C_C'][0], p['C_C_C'][0], i), lac(p['C_C_C'][1], p['C_C_C'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), lac(p['H_C_H'][0], p['Close'][0], i), lac(p['H_C_H'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), lac(p['C_C_H'][0], p['Dritt'][0], i), lac(p['C_C_H'][1], p['Dritt'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd1), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', metcar, methyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), lac(p['C_C_C_H'][0], p['0_1'][0], i), lac(p['C_C_C_H'][1], p['0_1'][1], i), lac(p['C_C_C_H'][2], p['0_1'][2], i), lac(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][2], p['0_C'][2], i), lac(p['CT'][3], p['0_C'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, lac(p['0_H'][2], p['HC'][2], i), lac(p['0_H'][3], p['HC'][3], i))

		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', alcoxy), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', metcar), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', hydhyd2), lac(p['H_C_H'][0], p['H_C_H'][0], i), lac(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(metcar, 'CT', alcoxy), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', metcar, methyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', alcoxy, alchyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', metcar), lac(p['C_C_O_H_2'][0], p['0_3'][0], i), lac(p['C_C_O_H_2'][1], p['0_3'][1], i), lac(p['C_C_O_H_2'][2], p['0_3'][2], i), lac(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', metcar), lac(p['C_C_O_H_1'][0], p['0_2'][0], i), lac(p['C_C_O_H_1'][1], p['0_2'][1], i), lac(p['C_C_O_H_1'][2], p['0_2'][2], i), lac(p['C_C_O_H_1'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', metcar, methyd), lac(p['C_C_O_H_2'][0], p['0_3'][0], i), lac(p['C_C_O_H_2'][1], p['0_3'][1], i), lac(p['C_C_O_H_2'][2], p['0_3'][2], i), lac(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', metcar, methyd), lac(p['0_2'][0], p['0_2'][0], i), lac(p['0_2'][1], p['0_2'][1], i), lac(p['0_2'][2], p['0_2'][2], i), lac(p['0_2'][3], p['0_2'][3], i))

		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), alcoxy, lac(p['OH'][0], p['0_O'][0], i), lac(p['OH'][1], p['0_O'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), alchyd, lac(p['HO'][0], p['0_H'][0], i), lac(p['HO'][1], p['0_H'][1], i))
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, lac(p['0_H'][0], p['HC'][0], i), lac(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', alcoxy), lac(p['CT_OH'][0], p['OH_mH'][0], i), lac(p['CT_OH'][1], p['OH_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd2), lac(p['HC_sO'][0], p['CT_HC'][0], i), lac(p['HC_sO'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(alcoxy, alchyd), lac(p['OH_HO'][0], p['HO_mH'][0], i), lac(p['OH_HO'][1], p['HO_mH'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', alcoxy), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', alcoxy), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd2), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', alcoxy, alchyd), lac(p['C_O_H'][0], p['Dritt'][0], i), lac(p['C_O_H'][1], p['Dritt'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', alcoxy, alchyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', 'CT', 'H1'), lac(p['C_C_O_H_2'][0], p['0_3'][0], i), lac(p['C_C_O_H_2'][1], p['0_3'][1], i), lac(p['C_C_O_H_2'][2], p['0_3'][2], i), lac(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', 'CT', 'H1'), lac(p['C_C_O_H_1'][0], p['0_2'][0], i), lac(p['C_C_O_H_1'][1], p['0_2'][1], i), lac(p['C_C_O_H_1'][2], p['0_2'][2], i), lac(p['C_C_O_H_1'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', 'CT'), lac(p['C_C_O_H_2'][0], p['0_3'][0], i), lac(p['C_C_O_H_2'][1], p['0_3'][1], i), lac(p['C_C_O_H_2'][2], p['0_3'][2], i), lac(p['C_C_O_H_2'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alchyd, alcoxy, 'CT', 'CT'), lac(p['C_C_O_H_1'][0], p['0_2'][0], i), lac(p['C_C_O_H_1'][1], p['0_2'][1], i), lac(p['C_C_O_H_1'][2], p['0_2'][2], i), lac(p['C_C_O_H_1'][3], p['0_2'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), alcoxy, lac(p['OH'][2], p['0_O'][2], i), lac(p['OH'][3], p['0_O'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), alchyd, lac(p['HO'][2], p['0_H'][2], i), lac(p['HO'][3], p['0_H'][3], i))
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, lac(p['0_H'][2], p['HC'][2], i), lac(p['0_H'][3], p['HC'][3], i))

		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', hydhyd1), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', hydhyd2), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', metcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', alcoxy), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amicar, 'CT', hydhyd1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amicar, 'CT', hydhyd2), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amicar, 'CT', metcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amicar, 'CT', alcoxy), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd3, 'CT', metcar, methyd), cal(p['0_1'][0], p['X_C_C_X'][0], i), cal(p['0_1'][1], p['X_C_C_X'][1], i), cal(p['0_1'][2], p['X_C_C_X'][2], i), cal(p['0_1'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd3, 'CT', alcoxy, alchyd), cal(p['0_5'][0], p['X_C_O_X'][0], i), cal(p['0_5'][1], p['X_C_O_X'][1], i), cal(p['0_5'][2], p['X_C_O_X'][2], i), cal(p['0_5'][3], p['X_C_O_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', amicar, aminit), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', amicar, amioxy), cal(p['H_C_C_O_1'][0], p['0_10'][0], i), cal(p['H_C_C_O_1'][1], p['0_10'][1], i), cal(p['H_C_C_O_1'][2], p['0_10'][2], i), cal(p['H_C_C_O_1'][3], p['0_10'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', amicar, amioxy), cal(p['H_C_C_O_2'][0], p['0_8'][0], i), cal(p['H_C_C_O_2'][1], p['0_8'][1], i), cal(p['H_C_C_O_2'][2], p['0_8'][2], i), cal(p['H_C_C_O_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', amicar, amioxy), cal(p['H_C_C_O_3'][0], p['0_9'][0], i), cal(p['H_C_C_O_3'][1], p['0_9'][1], i), cal(p['H_C_C_O_3'][2], p['0_9'][2], i), cal(p['H_C_C_O_3'][3], p['0_9'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', amicar, aminit), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', amicar, amioxy), cal(p['H_C_C_O_1'][0], p['0_10'][0], i), cal(p['H_C_C_O_1'][1], p['0_10'][1], i), cal(p['H_C_C_O_1'][2], p['0_10'][2], i), cal(p['H_C_C_O_1'][3], p['0_10'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', amicar, amioxy), cal(p['H_C_C_O_2'][0], p['0_8'][0], i), cal(p['H_C_C_O_2'][1], p['0_8'][1], i), cal(p['H_C_C_O_2'][2], p['0_8'][2], i), cal(p['H_C_C_O_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', amicar, amioxy), cal(p['H_C_C_O_3'][0], p['0_9'][0], i), cal(p['H_C_C_O_3'][1], p['0_9'][1], i), cal(p['H_C_C_O_3'][2], p['0_9'][2], i), cal(p['H_C_C_O_3'][3], p['0_9'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar, 'CT', amicar, amioxy), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar, 'CT', amicar, aminit), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', amicar, amioxy), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(alcoxy, 'CT', amicar, aminit), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amicar, 'CT', metcar, methyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amicar, 'CT', alcoxy, alchyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))

		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), amicar, cal(p['C'][0], p['0_C'][0], i), cal(p['C'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), amioxy, cal(p['O2'][0], p['0_O'][0], i), cal(p['O2'][1], p['0_O'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), aminit, cal(p['NA'][0], p['0_N'][0], i), cal(p['NA'][1], p['0_N'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), amihyd, cal(p['H'][0], p['0_H'][0], i), cal(p['H'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd3, cal(p['0_H'][0], p['H1'][0], i), cal(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', amicar), cal(p['CT_C'][0], p['CT_mH'][0], i), cal(p['CT_C'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd3), cal(p['HC_sC2'][0], p['CT_HC'][0], i), cal(p['HC_sC2'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(amicar, amioxy), cal(p['C_O'][0], p['O_mH'][0], i), cal(p['C_O'][1], p['O_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(amicar, aminit), cal(p['C_N'][0], p['N_mH'][0], i), cal(p['C_N'][1], p['N_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(aminit, amihyd), cal(p['NA_H'][0], p['H_mHC'][0], i), cal(p['NA_H'][1], p['H_mHC'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', amicar, amioxy), cal(p['C_C_O'][0], p['Dritt'][0], i), cal(p['C_C_O'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', amicar, aminit), cal(p['C_C_N'][0], p['Dritt'][0], i), cal(p['C_C_N'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amioxy, amicar, aminit), cal(p['O_C_N'][0], p['Close'][0], i), cal(p['O_C_N'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', amicar), cal(p['CT_CT_C'][0], p['C_C_H'][0], i), cal(p['CT_CT_C'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', amicar), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd3), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amicar, aminit, amihyd), cal(p['F_CA_CA_HA'][0], p['Dritt'][0], i), cal(p['F_CA_CA_HA'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(amihyd, aminit, amihyd), cal(p['H_N_H'][0], p['Close'][0], i), cal(p['H_N_H'][1], p['Close'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, amioxy), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), cal(p['C_C_C_N_1'][0], p['0_3'][0], i), cal(p['C_C_C_N_1'][1], p['0_3'][1], i), cal(p['C_C_C_N_1'][2], p['0_3'][2], i), cal(p['C_C_C_N_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), cal(p['C_C_C_N_2'][0], p['0_8'][0], i), cal(p['C_C_C_N_2'][1], p['0_8'][1], i), cal(p['C_C_C_N_2'][2], p['0_8'][2], i), cal(p['C_C_C_N_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), cal(p['C_C_C_N_3'][0], p['0_2'][0], i), cal(p['C_C_C_N_3'][1], p['0_2'][1], i), cal(p['C_C_C_N_3'][2], p['0_2'][2], i), cal(p['C_C_C_N_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', amicar, aminit), cal(p['C_C_C_N_4'][0], p['0_7'][0], i), cal(p['C_C_C_N_4'][1], p['0_7'][1], i), cal(p['C_C_C_N_4'][2], p['0_7'][2], i), cal(p['C_C_C_N_4'][3], p['0_7'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd3, 'CT', amicar, amioxy), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd3, 'CT', amicar, aminit), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amihyd, aminit, amicar, amioxy), cal(p['O_C_N_H_1'][0], p['0_3'][0], i), cal(p['O_C_N_H_1'][1], p['0_3'][1], i), cal(p['O_C_N_H_1'][2], p['0_3'][2], i), cal(p['O_C_N_H_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amihyd, aminit, amicar, amioxy), cal(p['O_C_N_H_2'][0], p['0_11'][0], i), cal(p['O_C_N_H_2'][1], p['0_11'][1], i), cal(p['O_C_N_H_2'][2], p['0_11'][2], i), cal(p['O_C_N_H_2'][3], p['0_11'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amihyd, aminit, amicar, 'CT'), cal(p['X_C_N_X'][0], p['Ring_0'][0], i), cal(p['X_C_N_X'][1], p['Ring_0'][1], i), cal(p['X_C_N_X'][2], p['Ring_0'][2], i), cal(p['X_C_N_X'][3], p['Ring_0'][3], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', amicar, amioxy), cal(p['Car_imp'][0], p['Imp_0'][0], i), cal(p['Car_imp'][1], p['Imp_0'][1], i), cal(p['Car_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', 'X ', aminit, amihyd), cal(p['Ami_imp'][0], p['Imp_0'][0], i), cal(p['Ami_imp'][1], p['Imp_0'][1], i), cal(p['Ami_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amicar, 'CT', aminit, amihyd), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(amicar, 'CT', aminit, amioxy), cal(p['Ring_imp'][0], p['Imp_0'][0], i), cal(p['Ring_imp'][1], p['Imp_0'][1], i), cal(p['Ring_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), amicar, cal(p['C'][2], p['0_C'][2], i), cal(p['C'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), amioxy, cal(p['O2'][2], p['0_O'][2], i), cal(p['O2'][3], p['0_O'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), aminit, cal(p['NA'][2], p['0_N'][2], i), cal(p['NA'][3], p['0_N'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), amihyd, cal(p['H'][2], p['0_H'][2], i), cal(p['H'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd3, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))