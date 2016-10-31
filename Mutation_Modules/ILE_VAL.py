# ILE to VAL Mutation

import Frcmod_creator
import PDBHandler
import Leapy
from ParmedTools.ParmedActions import *
from chemistry.amber.readparm import *

def parmed_command(vxi='VXI'):
	bc = {}
        with open('Param_files/AminoAcid/ILE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/VAL.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HG11 :{}@HD11 0 0'.format(vxi, vxi)).execute()
                change(parm, 'charge', ':{}@N'.format(vxi), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vxi), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vxi), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vxi), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vxi), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB'.format(vxi), bc['HB']+((fc['HB']-bc['HB'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG2'.format(vxi), bc['CG2']+((fc['CG2']-bc['CG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG21'.format(vxi), bc['HG21']+((fc['HG21']-bc['HG21'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG22'.format(vxi), bc['HG22']+((fc['HG22']-bc['HG22'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG23'.format(vxi), bc['HG23']+((fc['HG23']-bc['HG23'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG1'.format(vxi), bc['CG1']+((fc['CG1']-bc['CG1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG11'.format(vxi), (fc['HG11']/10)*i).execute()
                change(parm, 'charge', ':{}@HG12'.format(vxi), bc['HG12']+((fc['HG12']-bc['HG12'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG13'.format(vxi), bc['HG13']+((fc['HG13']-bc['HG13'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD1'.format(vxi), bc['CG2']-(bc['CG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HD11'.format(vxi), bc['HG21']-(bc['HG21']/10)*i).execute()
                change(parm, 'charge', ':{}@HD12'.format(vxi), bc['HG22']-(bc['HG22']/10)*i).execute()
                change(parm, 'charge', ':{}@HD13'.format(vxi), bc['HG23']-(bc['HG23']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vxi), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vxi), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi(struct, out, aa, vxi='VXI'):
        struct.residue_dict[aa].set_resname(vxi)
        CD1 = struct.residue_dict[aa].atom_dict['CD1']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'CG1' and res.get_resname() == vxi:
                        	pdb.write(atom.formatted())
                        	pdb.write(atom.superimposed1('HG11', CD1))
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

def lib_make(ff, outputfile, vxi='VXI', metcar='dc', methyd='dh', hydhyd='mh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source leaprc.%s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ILE-VAL.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "C"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "H"\n'%vxi)
	ctrl.write('set %s.1.19 element "C"\n'%vxi)
	ctrl.write('set %s.1.20 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB"\n'%vxi)
	ctrl.write('set %s.1.7 name "CG2"\n'%vxi)
	ctrl.write('set %s.1.8 name "HG21"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG22"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG23"\n'%vxi)
	ctrl.write('set %s.1.11 name "CG1"\n'%vxi)
	ctrl.write('set %s.1.12 name "HG11"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG12"\n'%vxi)
	ctrl.write('set %s.1.14 name "HG13"\n'%vxi)
	ctrl.write('set %s.1.15 name "CD1"\n'%vxi)
	ctrl.write('set %s.1.16 name "HD11"\n'%vxi)
	ctrl.write('set %s.1.17 name "HD12"\n'%vxi)
	ctrl.write('set %s.1.18 name "HD13"\n'%vxi)
	ctrl.write('set %s.1.19 name "C"\n'%vxi)
	ctrl.write('set %s.1.20 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "CT"\n'%vxi)
	ctrl.write('set %s.1.8 type "HC"\n'%vxi)
	ctrl.write('set %s.1.9 type "HC"\n'%vxi)
	ctrl.write('set %s.1.10 type "HC"\n'%vxi)
	ctrl.write('set %s.1.11 type "CT"\n'%vxi)
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, hydhyd))
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.16 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.17 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.19 type "C"\n'%vxi)
	ctrl.write('set %s.1.20 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.7 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.19 %s.1.20\n'%(vxi, vxi))
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

def stock_add_to_all(metcar='dc', methyd='dh', hydhyd='mh'):
	Frcmod_creator.make_hyb()
	Frcmod_creator.TYPE_insert(metcar, 'C', 'sp3')
	Frcmod_creator.TYPE_insert(methyd, 'H', 'sp3')
	Frcmod_creator.TYPE_insert(hydhyd, 'H', 'sp3')
	bp = {}
	with open('Param_files/Stock/MGrown.param', 'r') as b:
		data = b.readlines()[1:]
	for line in data:
		bp[line.split()[0]] = []
		for point in line.split()[1:]:
			bp[line.split()[0]].append(float(point))
	b.close()
	fp = {}
	with open('Param_files/Stock/MSmall.param', 'r') as b:
		data = b.readlines()[1:]
	for line in data:
		fp[line.split()[0]] = []
		for point in line.split()[1:]:
			fp[line.split()[0]].append(float(point))
	b.close()
	for i in range(11):
		a = i*10
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, bp['dc'][0]+((fp['dc'][0]-bp['dc'][0])/10)*i, bp['dc'][1]+((fp['dc'][1]-bp['dc'][1])/10)*i)
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, bp['dh'][0]+((fp['dh'][0]-bp['dh'][0])/10)*i, bp['dh'][1]+((fp['dh'][1]-bp['dh'][1])/10)*i)
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, bp['mh'][0]+((fp['mh'][0]-bp['mh'][0])/10)*i, bp['mh'][1]+((fp['mh'][1]-bp['mh'][1])/10)*i)
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), bp['CT-dc'][0]+((fp['CT-dc'][0]-bp['CT-dc'][0])/10)*i, bp['CT-dc'][1]+((fp['CT-dc'][1]-bp['CT-dc'][1])/10)*i)
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), bp['CT-mh'][0]+((fp['CT-mh'][0]-bp['CT-mh'][0])/10)*i, bp['CT-mh'][1]+((fp['CT-mh'][1]-bp['CT-mh'][1])/10)*i)
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), bp['dc-dh'][0]+((fp['dc-dh'][0]-bp['dc-dh'][0])/10)*i, bp['dc-dh'][1]+((fp['dc-dh'][1]-bp['dc-dh'][1])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), bp['CT-dc-dh'][0]+((fp['CT-dc-dh'][0]-bp['CT-dc-dh'][0])/10)*i, bp['CT-dc-dh'][1]+((fp['CT-dc-dh'][1]-bp['CT-dc-dh'][1])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), bp['dh-dc-dh'][0]+((fp['dh-dc-dh'][0]-bp['dh-dc-dh'][0])/10)*i, bp['dh-dc-dh'][1]+((fp['dh-dc-dh'][1]-bp['dh-dc-dh'][1])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', metcar), bp['HC-CT-dc'][0]+((fp['HC-CT-dc'][0]-bp['HC-CT-dc'][0])/10)*i, bp['HC-CT-dc'][1]+((fp['HC-CT-dc'][1]-bp['HC-CT-dc'][1])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', hydhyd), bp['HC-CT-mh'][0]+((fp['HC-CT-mh'][0]-bp['HC-CT-mh'][0])/10)*i, bp['HC-CT-mh'][1]+((fp['HC-CT-mh'][1]-bp['HC-CT-mh'][1])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', metcar), bp['mh-CT-dc'][0]+((fp['mh-CT-dc'][0]-bp['mh-CT-dc'][0])/10)*i, bp['mh-CT-dc'][1]+((fp['mh-CT-dc'][1]-bp['mh-CT-dc'][1])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), bp['CT-CT-dc'][0]+((fp['CT-CT-dc'][0]-bp['CT-CT-dc'][0])/10)*i, bp['CT-CT-dc'][1]+((fp['CT-CT-dc'][1]-bp['CT-CT-dc'][1])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd), bp['CT-CT-mh'][0]+((fp['CT-CT-mh'][0]-bp['CT-CT-mh'][0])/10)*i, bp['CT-CT-mh'][1]+((fp['CT-CT-mh'][1]-bp['CT-CT-mh'][1])/10)*i)
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), bp['CT-CT-dc-dh'][0]+((fp['CT-CT-dc-dh'][0]-bp['CT-CT-dc-dh'][0])/10)*i, bp['CT-CT-dc-dh'][1]+((fp['CT-CT-dc-dh'][1]-bp['CT-CT-dc-dh'][1])/10)*i, bp['CT-CT-dc-dh'][2]+((fp['CT-CT-dc-dh'][2]-bp['CT-CT-dc-dh'][2])/10)*i, bp['CT-CT-dc-dh'][3]+((fp['CT-CT-dc-dh'][3]-bp['CT-CT-dc-dh'][3])/10)*i)
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', metcar, methyd), bp['HC-CT-dc-dh'][0]+((fp['HC-CT-dc-dh'][0]-bp['HC-CT-dc-dh'][0])/10)*i, bp['CT-CT-dc-dh'][1]+((fp['HC-CT-dc-dh'][1]-bp['HC-CT-dc-dh'][1])/10)*i, bp['HC-CT-dc-dh'][2]+((fp['HC-CT-dc-dh'][2]-bp['HC-CT-dc-dh'][2])/10)*i, bp['HC-CT-dc-dh'][3]+((fp['HC-CT-dc-dh'][3]-bp['HC-CT-dc-dh'][3])/10)*i)
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', metcar, methyd), bp['mh-CT-dc-dh'][0]+((fp['mh-CT-dc-dh'][0]-bp['mh-CT-dc-dh'][0])/10)*i, bp['CT-CT-dc-dh'][1]+((fp['mh-CT-dc-dh'][1]-bp['mh-CT-dc-dh'][1])/10)*i, bp['mh-CT-dc-dh'][2]+((fp['mh-CT-dc-dh'][2]-bp['mh-CT-dc-dh'][2])/10)*i, bp['mh-CT-dc-dh'][3]+((fp['mh-CT-dc-dh'][3]-bp['mh-CT-dc-dh'][3])/10)*i)
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, bp['dc'][2]+((fp['dc'][2]-bp['dc'][2])/10)*i, bp['dc'][3]+((fp['dc'][3]-bp['dc'][3])/10)*i)
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, bp['dh'][2]+((fp['dh'][2]-bp['dh'][2])/10)*i, bp['dh'][3]+((fp['dh'][3]-bp['dh'][3])/10)*i)
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, bp['mh'][2]+((fp['mh'][2]-bp['mh'][2])/10)*i, bp['mh'][3]+((fp['mh'][3]-bp['mh'][3])/10)*i)
