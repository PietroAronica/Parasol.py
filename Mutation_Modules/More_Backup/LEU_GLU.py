# LEU to GLU Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
	bc = {}
        with open('Param_files/AminoAcid/LEU.param', 'r') as b:
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
                changeLJPair(parm, ':{}@HG :{}@OE1 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HG2 :{}@HD11 0 0'.format(vxi, vxi)).execute()
                changeLJPair(parm, ':{}@HG3 :{}@HD21 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vxi), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vxi), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vxi), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG'.format(vxi), bc['HG']-(bc['HG']/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vxi), (fc['HG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vxi), (fc['HG3']/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vxi), (fc['CD']/10)*i*i/10).execute()
                change(parm, 'charge', ':{}@OE1'.format(vxi), (fc['OE1']/10)*i*i/10).execute()
                change(parm, 'charge', ':{}@OE2'.format(vxi), (fc['OE2']/10)*i*i/10).execute()
                change(parm, 'charge', ':{}@CD1'.format(vxi), bc['CD1']-(bc['HD11']/10)*i).execute()
                change(parm, 'charge', ':{}@HD11'.format(vxi), bc['HD11']-(bc['HD11']/10)*i).execute()
                change(parm, 'charge', ':{}@HD12'.format(vxi), bc['HD12']-(bc['HD12']/10)*i).execute()
                change(parm, 'charge', ':{}@HD13'.format(vxi), bc['HD13']-(bc['HD13']/10)*i).execute()
                change(parm, 'charge', ':{}@CD2'.format(vxi), bc['CD2']-(bc['HD21']/10)*i).execute()
                change(parm, 'charge', ':{}@HD21'.format(vxi), bc['HD21']-(bc['HD21']/10)*i).execute()
                change(parm, 'charge', ':{}@HD22'.format(vxi), bc['HD22']-(bc['HD22']/10)*i).execute()
                change(parm, 'charge', ':{}@HD23'.format(vxi), bc['HD23']-(bc['HD23']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
                d = netCharge(parm).execute()
                change(parm, 'charge', ':PN', '{:.3f}'.format(-d)).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CG = struct.residue_dict[aa].atom_dict['CG']
        HG = struct.residue_dict[aa].atom_dict['HG']
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
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.superimposed1('HG2', CD1))
                        	pdb.write(atom.superimposed1('HG3', CD2))
                        	pdb.write(atom.halfway_between('CD', CG, HG))
                                pdb.write(atom.superimposed1('OE1', HG))
                                pdb.write(atom.superimposed2('OE2', HG))
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

def lib_make(ff, outputfile, vxi='VXI', metcar1='1c', methyd1='1h', hydhyd1='xh', metcar2='2c', methyd2='2h', hydhyd2='yh', carcar='cc', caroxy='co', hydhyd3='sh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/LEU-GLU.pdb\n"%vxi)
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
	ctrl.write('set %s.1.13 element "O"\n'%vxi)
	ctrl.write('set %s.1.14 element "O"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
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
	ctrl.write('set %s.1.9 name "HG"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.12 name "CD"\n'%vxi)
	ctrl.write('set %s.1.13 name "OE1"\n'%vxi)
	ctrl.write('set %s.1.14 name "OE2"\n'%vxi)
	ctrl.write('set %s.1.15 name "CD1"\n'%vxi)
	ctrl.write('set %s.1.16 name "HD11"\n'%vxi)
	ctrl.write('set %s.1.17 name "HD12"\n'%vxi)
	ctrl.write('set %s.1.18 name "HD13"\n'%vxi)
	ctrl.write('set %s.1.19 name "CD2"\n'%vxi)
	ctrl.write('set %s.1.20 name "HD21"\n'%vxi)
	ctrl.write('set %s.1.21 name "HD22"\n'%vxi)
	ctrl.write('set %s.1.22 name "HD23"\n'%vxi)
	ctrl.write('set %s.1.23 name "C"\n'%vxi)
	ctrl.write('set %s.1.24 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%(vxi))
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "CT"\n'%vxi)
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, hydhyd3))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, hydhyd1))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, hydhyd2))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, carcar))
	ctrl.write('set %s.1.13 type "%s"\n'%(vxi, caroxy))
	ctrl.write('set %s.1.14 type "%s"\n'%(vxi, caroxy))
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, metcar2))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, methyd2))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, methyd2))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, methyd2))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, metcar1))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, methyd1))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, methyd1))
	ctrl.write('set %s.1.22 type "%s"\n'%(vxi, methyd1))
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
	ctrl.write('bond %s.1.8 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.22\n'%(vxi, vxi))
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

def cal2(x, y, i):
        num = y+((x-y)/10)*i
        return num

def stock_add_to_all(vxi='VXI', metcar1='1c', methyd1='1h', hydhyd1='xh', metcar2='2c', methyd2='2h', hydhyd2='yh', carcar='cc', caroxy='co', hydhyd3='sh'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(metcar1, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(methyd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd1, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(metcar2, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(methyd2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd2, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(carcar, 'C', 'sp2')
	Frcmod_creator.TYPE_insert(caroxy, 'O', 'sp2')
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
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar1, cal(p['CT'][0], p['0_C'][0], i), cal(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd1, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar1), cal(p['CT_CT'][0], p['CT_mH'][0], i), cal(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd1), cal(p['HC_sC'][0], p['CT_HC'][0], i), cal(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar1, methyd1), cal(p['CT_HC'][0], p['HC_mH'][0], i), cal(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar1, methyd1), cal(p['C_C_H'][0], p['Dritt'][0], i), cal(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd1, metcar1, methyd1), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(metcar2, 'CT', metcar1), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', metcar1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', metcar2), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', hydhyd2), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar1), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carcar, 'CT', metcar1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carcar, 'CT', hydhyd1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd1, 'CT', metcar1), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar1, methyd1), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', metcar1, methyd1), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', metcar1, methyd1), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar1, cal(p['CT'][2], p['0_C'][2], i), cal(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd1, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd1, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar2, cal(p['CT'][0], p['0_C'][0], i), cal(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd2, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, cal(p['0_H'][0], p['HC'][0], i), cal(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar2), cal(p['CT_CT'][0], p['CT_mH'][0], i), cal(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd2), cal(p['HC_sC'][0], p['CT_HC'][0], i), cal(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar2, methyd2), cal(p['CT_HC'][0], p['HC_mH'][0], i), cal(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar2, methyd2), cal(p['C_C_H'][0], p['Dritt'][0], i), cal(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd2, metcar2, methyd2), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar2), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', metcar1), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', metcar2), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carcar, 'CT', metcar2), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', hydhyd2), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', hydhyd1), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(carcar, 'CT', hydhyd2), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd2, 'CT', metcar2), cal(p['Close'][0], p['Close'][0], i), cal(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd2), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar2, methyd2), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar1, 'CT', metcar2, methyd2), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar2, 'CT', metcar1, methyd1), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd3, 'CT', metcar2, methyd2), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd3, 'CT', metcar1, methyd1), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', metcar2, methyd2), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', metcar1, methyd1), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', metcar2, methyd2), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carcar, 'CT', metcar2, methyd2), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(carcar, 'CT', metcar1, methyd1), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar2, cal(p['CT'][2], p['0_C'][2], i), cal(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd2, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd2, cal(p['0_H'][2], p['HC'][2], i), cal(p['0_H'][3], p['HC'][3], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), carcar, cal2(p['C'][0], p['0_C'][0], i), cal2(p['C'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), caroxy, cal2(p['O2'][0], p['0_O'][0], i), cal2(p['O2'][1], p['0_O'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd3, cal2(p['0_H'][0], p['HC'][0], i), cal2(p['0_H'][1], p['HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', carcar), cal2(p['CT_C'][0], p['CT_mH'][0], i), cal2(p['CT_C'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd3), cal2(p['HC_sC2'][0], p['CT_HC'][0], i), cal2(p['HC_sC2'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(carcar, caroxy), cal2(p['C_O2'][0], p['O2_mH'][0], i), cal2(p['C_O2'][1], p['O2_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', carcar, caroxy), cal2(p['C_C_O2'][0], p['Dritt'][0], i), cal2(p['C_C_O2'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(caroxy, carcar, caroxy), cal2(p['O2_C_O2'][0], p['Close'][0], i), cal2(p['O2_C_O2'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', carcar), cal2(p['CT_CT_C'][0], p['C_C_H'][0], i), cal2(p['CT_CT_C'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', carcar), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', hydhyd3), cal2(p['H_C_H'][0], p['H_C_H'][0], i), cal2(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd3, 'CT', carcar), cal2(p['Close'][0], p['Close'][0], i), cal2(p['Close'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd3), cal2(p['C_C_H'][0], p['C_C_H'][0], i), cal2(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', carcar, caroxy), cal2(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal2(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal2(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal2(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', carcar, caroxy), cal2(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal2(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal2(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal2(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd2, 'CT', carcar, caroxy), cal2(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal2(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal2(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal2(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd1, 'CT', carcar, caroxy), cal2(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal2(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal2(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal2(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd3, 'CT', carcar, caroxy), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar2, 'CT', carcar, caroxy), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(metcar1, 'CT', carcar, caroxy), cal2(p['0_Dihe'][0], p['0_Dihe'][0], i), cal2(p['0_Dihe'][1], p['0_Dihe'][1], i), cal2(p['0_Dihe'][2], p['0_Dihe'][2], i), cal2(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('X ', caroxy, carcar, caroxy), cal2(p['Car_imp'][0], p['Imp_0'][0], i), cal2(p['Car_imp'][1], p['Imp_0'][1], i), cal2(p['Car_imp'][2], p['Imp_0'][2], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), carcar, cal2(p['C'][2], p['0_C'][2], i), cal2(p['C'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), caroxy, cal2(p['O2'][2], p['0_O'][2], i), cal2(p['O2'][3], p['0_O'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd3, cal2(p['0_H'][2], p['HC'][2], i), cal2(p['0_H'][3], p['HC'][3], i))
