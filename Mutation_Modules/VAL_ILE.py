# VAL to ILE Mutation

import Frcmod_creator
import PDBHandler

def makevxi(struct, out, vxi):
        struct.residue_dict[vxi].set_resname('VXI')
        CG2 = struct.residue_dict[vxi].atom_dict['CG2']
        HG21 = struct.residue_dict[vxi].atom_dict['HG21']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'CG1' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('CG2'))
			elif atom.get_name() == 'HG11' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('HG21'))
			elif atom.get_name() == 'HG12' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('HG22'))
			elif atom.get_name() == 'HG13' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('HG23'))
			elif atom.get_name() == 'CG2' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('CG1'))
			elif atom.get_name() == 'HG21' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('HG11'))
                        	pdb.write(atom.halfway_between('CD1', CG2, HG21))
                        	pdb.write(atom.superimposed1('HD11', HG21))
                        	pdb.write(atom.superimposed2('HD12', HG21))
                        	pdb.write(atom.superimposed3('HD13', HG21))
			elif atom.get_name() == 'HG22' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('HG12'))
			elif atom.get_name() == 'HG23' and res.get_resname() == 'VXI':
                        	pdb.write(atom.change_name('HG13'))
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


def all_make():
	for i in range(0,110,10):
		Frcmod_creator.make ('{}_{}.frcmod'.format(i, 100-i))

def stock_add_to_all():
	begparam = {}
	with open('Param_files/Stock/Methyl_Small.param', 'r') as b:
		data = b.readlines()[1:]
	for line in data:
		key, value = line.split("=")
		begparam[key.strip()] = float(value.strip())
	b.close()
	finparam = {}
	with open('Param_files/Stock/Methyl_Grown.param', 'r') as b:
		data = b.readlines()[1:]
	for line in data:
		key, value = line.split("=")
		finparam[key.strip()] = float(value.strip())
	b.close()
	for i in range(11):
		a = i*10
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), 'dc', begparam['dc_mass']+((finparam['dc_mass']-begparam['dc_mass'])/10)*i, begparam['dc_pola']+((finparam['dc_pola']-begparam['dc_pola'])/10)*i)
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), 'dh', begparam['dh_mass']+((finparam['dh_mass']-begparam['dh_mass'])/10)*i, begparam['dh_pola']+((finparam['dh_pola']-begparam['dh_pola'])/10)*i)
		Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), 'mh', begparam['mh_mass']+((finparam['mh_mass']-begparam['mh_mass'])/10)*i, begparam['mh_pola']+((finparam['mh_pola']-begparam['mh_pola'])/10)*i)
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), 'CT-dc', begparam['CT-dc_fconst']+((finparam['CT-dc_fconst']-begparam['CT-dc_fconst'])/10)*i, begparam['CT-dc_length']+((finparam['CT-dc_length']-begparam['CT-dc_length'])/10)*i)
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), 'dc-dh', begparam['dc-dh_fconst']+((finparam['dc-dh_fconst']-begparam['dc-dh_fconst'])/10)*i, begparam['dc-dh_length']+((finparam['dc-dh_length']-begparam['dc-dh_length'])/10)*i)
		Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), 'CT-mh', begparam['CT-mh_fconst']+((finparam['CT-mh_fconst']-begparam['CT-mh_fconst'])/10)*i, begparam['CT-mh_length']+((finparam['CT-mh_length']-begparam['CT-mh_length'])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), 'CT-dc-dh', begparam['CT-dc-dh_fcnst']+((finparam['CT-dc-dh_fcnst']-begparam['CT-dc-dh_fcnst'])/10)*i, begparam['CT-dc-dh_width']+((finparam['CT-dc-dh_width']-begparam['CT-dc-dh_width'])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), 'dh-dc-dh', begparam['dh-dc-dh_fcnst']+((finparam['dh-dc-dh_fcnst']-begparam['dh-dc-dh_fcnst'])/10)*i, begparam['dh-dc-dh_width']+((finparam['dh-dc-dh_width']-begparam['dh-dc-dh_width'])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), 'HC-CT-dc', begparam['HC-CT-dc_fcnst']+((finparam['HC-CT-dc_fcnst']-begparam['HC-CT-dc_fcnst'])/10)*i, begparam['HC-CT-dc_width']+((finparam['HC-CT-dc_width']-begparam['HC-CT-dc_width'])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), 'HC-CT-mh', begparam['HC-CT-mh_fcnst']+((finparam['HC-CT-mh_fcnst']-begparam['HC-CT-mh_fcnst'])/10)*i, begparam['HC-CT-mh_width']+((finparam['HC-CT-mh_width']-begparam['HC-CT-mh_width'])/10)*i)
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), 'CT-dc-dh', begparam['CT-dc-dh_fcnst']+((finparam['CT-dc-dh_fcnst']-begparam['CT-dc-dh_fcnst'])/10)*i, begparam['CT-dc-dh_width']+((finparam['CT-dc-dh_width']-begparam['CT-dc-dh_width'])/10)*i)
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), 'mh-CT-dc-dh', begparam['mh-CT-dc-dh_div']+((finparam['mh-CT-dc-dh_div']-begparam['mh-CT-dc-dh_div'])/10)*i, begparam['mh-CT-dc-dh_fcs']+((finparam['mh-CT-dc-dh_fcs']-begparam['mh-CT-dc-dh_fcs'])/10)*i, begparam['mh-CT-dc-dh_ang']+((finparam['mh-CT-dc-dh_ang']-begparam['mh-CT-dc-dh_ang'])/10)*i, begparam['mh-CT-dc-dh_per']+((finparam['mh-CT-dc-dh_per']-begparam['mh-CT-dc-dh_per'])/10)*i)
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), 'dc', begparam['dc_vdwr']+((finparam['dc_vdwr']-begparam['dc_vdwr'])/10)*i, begparam['dc_ptwd']+((finparam['dc_ptwd']-begparam['dc_ptwd'])/10)*i)
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), 'dh', begparam['dh_vdwr']+((finparam['dh_vdwr']-begparam['dh_vdwr'])/10)*i, begparam['dh_ptwd']+((finparam['dh_ptwd']-begparam['dh_ptwd'])/10)*i)
		Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), 'mh', begparam['mh_vdwr']+((finparam['mh_vdwr']-begparam['mh_vdwr'])/10)*i, begparam['mh_ptwd']+((finparam['mh_ptwd']-begparam['mh_ptwd'])/10)*i)
