import numpy
import Leapy
import os, sys
import math

_CRYST1_FORMAT_STRING = '{:6}{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}\n' 
_CONECT_FORMAT_STRING = '{:6}{:5d}{:5d}\n' 
_ATOM_FORMAT_STRING = '{:4} {:6d} {:4} {:3} {:5}     {:7.3f} {:7.3f} {:7.3f} {:5.2f} {:5.2f}\n' 

class Atom(object):
	def __init__(self, number, name, resname, resnumber, chainid, coord, occupancy, bfactor, element=None):
		self.number = int(number)
		self.name = name
		self.resname = resname
		self.resnumber = int(resnumber)
		self.chainid = chainid
		self.coord = coord
		self.x = self.coord[0]
		self.y = self.coord[1]
		self.z = self.coord[2]
		self.occupancy = occupancy
		self.bfactor = bfactor
		self.id = '{} of {}{}'.format(name, resname, resnumber)

	def __repr__(self):
		return self.name

	def get_name(self):
		return self.name

	def set_resname(self, resname):
		self.resname = resname

	def set_name(self, name):
		self.name = name

	def get_coord(self):
		return self.coord

	def calc_distance(self, atom):
		ax = self.get_coord()[0]
		ay = self.get_coord()[1]
		az = self.get_coord()[2]
		bx = atom.get_coord()[0]
		by = atom.get_coord()[1]
		bz = atom.get_coord()[2]
		distance = math.sqrt((bx - ax) ** 2 + (by - ay) ** 2 + (bz - az) ** 2)
		return distance

	def get_number(self):
		return self.number

	def get_id(self):
		return self.id

	def formatted(self):
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, self.name, self.resname, self.resnumber, self.x, self.y, self.z, self.occupancy, self.bfactor)

	def change_name(self, name):
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, self.x, self.y, self.z, self.occupancy, self.bfactor)

	def superimposed1(self, name, atom):
		x = atom.get_coord()[0] - 0.001
		y = atom.get_coord()[1] - 0.001
		z = atom.get_coord()[2] - 0.001
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def superimposed2(self, name, atom):
		x = atom.get_coord()[0] - 0.002
		y = atom.get_coord()[1] - 0.002
		z = atom.get_coord()[2] - 0.002
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def superimposed3(self, name, atom):
		x = atom.get_coord()[0] + 0.001
		y = atom.get_coord()[1] + 0.001
		z = atom.get_coord()[2] + 0.001
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def superimposed4(self, name, atom):
		x = atom.get_coord()[0] + 0.002
		y = atom.get_coord()[1] + 0.002
		z = atom.get_coord()[2] + 0.002
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def halfway_between(self, name, atom1, atom2):
		x = (atom2.get_coord()[0] + atom1.get_coord()[0])/2
		y = (atom2.get_coord()[1] + atom1.get_coord()[1])/2
		z = (atom2.get_coord()[2] + atom1.get_coord()[2])/2
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def halfway_between1(self, name, atom1, atom2):
		x = ((atom2.get_coord()[0] + atom1.get_coord()[0])/2) - 0.001
		y = ((atom2.get_coord()[1] + atom1.get_coord()[1])/2) - 0.001
		z = ((atom2.get_coord()[2] + atom1.get_coord()[2])/2) - 0.001
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def halfway_between2(self, name, atom1, atom2):
		x = ((atom2.get_coord()[0] + atom1.get_coord()[0])/2) + 0.001
		y = ((atom2.get_coord()[1] + atom1.get_coord()[1])/2) + 0.001
		z = ((atom2.get_coord()[2] + atom1.get_coord()[2])/2) + 0.001
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def halfway_after1(self, name, atom1, atom2):
		x = ((3*atom2.get_coord()[0] - atom1.get_coord()[0])/2) - 0.001
		y = ((3*atom2.get_coord()[1] - atom1.get_coord()[1])/2) - 0.001
		z = ((3*atom2.get_coord()[2] - atom1.get_coord()[2])/2) - 0.001
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def halfway_after2(self, name, atom1, atom2):
		x = ((3*atom2.get_coord()[0] - atom1.get_coord()[0])/2) + 0.001
		y = ((3*atom2.get_coord()[1] - atom1.get_coord()[1])/2) + 0.001
		z = ((3*atom2.get_coord()[2] - atom1.get_coord()[2])/2) + 0.001
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def thirdone_between(self, name, atom1, atom2):
		x = (atom2.get_coord()[0] + 2*atom1.get_coord()[0])/3
		y = (atom2.get_coord()[1] + 2*atom1.get_coord()[1])/3
		z = (atom2.get_coord()[2] + 2*atom1.get_coord()[2])/3
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

	def thirdtwo_between(self, name, atom1, atom2):
		x = (2*atom2.get_coord()[0] + atom1.get_coord()[0])/3
		y = (2*atom2.get_coord()[1] + atom1.get_coord()[1])/3
		z = (2*atom2.get_coord()[2] + atom1.get_coord()[2])/3
		return _ATOM_FORMAT_STRING.format('ATOM', self.number, name, self.resname, self.resnumber, x, y, z, self.occupancy, self.bfactor)

class Residue(object):
	def __init__(self, resname, resnumber):
		self.resname = resname
		self.resnumber = resnumber
		self.id = '{}{}'.format(resname, resnumber)
		self.atom_list = []
		self.atom_dict = {}

	def __repr__(self):
		return self.id

	def get_resname(self):
		return self.resname

	def set_resname(self, resname):
		self.resname = resname
		for atom in self.atom_list:
			atom.set_resname(resname)

	def get_resnumber(self):
		return self.resnumber

	def get_id(self):
		return self.id

	def get_atomlist(self):
		for t in self.atom_dict:
			print t

	def add_atom(self, atom):
		atom_name = atom.get_name()
		self.atom_list.append(atom)
		self.atom_dict[atom_name] = atom

class Structure(object):
	def __init__(self, name):
		self.name = name
		self.residue_list = []
		self.residue_dict = {}
		self.other_list = []
		self.other_dict = {}
		self.id = 'PDB read as {}'.format(name)

	def __repr__(self):
		return self.id

	def get_reslist(self):
		for t in self.residue_dict:
			print t
	
	def add_residue(self, residue):
		residue_number = residue.get_resnumber()
		self.residue_list.append(residue)
		self.residue_dict[residue_number] = residue

	def add_ter(self, ter):
		ter_position = ter.ter_position()
		self.other_list.append(ter)
		self.other_dict[ter_position] = ter

	def add_cryst1(self, cryst1):
		self.other_list.append(cryst1)
		self.other_dict['Cryst1'] = cryst1

	def add_conect(self, conect):
		conect_name = conect.__repr__()
		self.other_list.append(conect)
		self.other_dict[conect_name] = conect

class Ter_record(object):
	def __init__(self, position):
		self.position = position
		self.tag = '{}\n'.format('TER')
		self.id = 'Ter following residue {}'.format(position)

	def __repr__(self):
		return self.id

	def ter(self):
		return self.tag

	def ter_position(self):
		return self.position

class Cryst1_record(object):
	def __init__(self, xbox, ybox, zbox, xangle, yangle, zangle):
		self.xbox = xbox
		self.ybox = ybox
		self.zbox = zbox
		self.boxdim = numpy.array((xbox, ybox, zbox), "f")
		self.boxangles = numpy.array((xangle, yangle, zangle), "f")
		self.xangle = xangle
		self.yangle = yangle
		self.zangle = zangle
		self.id = 'Cryst1 Record'

	def __repr__(self):
		return self.id

	def formatted(self):
		return _CRYST1_FORMAT_STRING.format('CRYST1', self.xbox, self.ybox, self.zbox, self.xangle, self.yangle, self.zangle)

	def box_side(self):
		return self.xbox

	def box_angle(self):
		return self.xangle

	def box_dimensions(self):
		return self.boxdim

	def box_angles(self):
		return self.boxangles

class Conect_record(object):
	def __init__(self, atom1, atom2):
		self.atom1 = int(atom1)
		self.atom2 = int(atom2)
		self.id = 'Conect Record between atoms {} and {}'.format(atom1, atom2)

	def __repr__(self):
		return self.id

	def formatted(self):
		return _CONECT_FORMAT_STRING.format('CONECT', self.atom1, self.atom2)

	def atoms(self):
		return self.atom1, self.atom2

def readpdb(file):
	structure = Structure('curr')
	pdb = open(file, 'r')
	data = pdb.readlines()
	current_resnumber = None
	for line in data:
		line = line.rstrip('\n')
		record_type = line [0:6].strip()
		if record_type == "ATOM" or record_type == "HETATM":
			number = line[6:11].strip()
			name = line[12:16].strip()
			resname = line[17:20].strip()
			chainid = line[21].strip()
			resnumber = int(line[22:27].split()[0])
			x = float(line[30:38])
			y = float(line[38:46])
			z = float(line[46:54])
			coord = numpy.array((x, y, z), "f")
			try:
				occupancy = float(line[54:60])
			except:
				occupancy = 1.00
			try:
				bfactor = float(line[60:66])
			except:
				bfactor = 0.00
			curr_atom = Atom(number, name, resname, resnumber, chainid, coord, occupancy, bfactor)
			try:
				if curr_residue.get_resnumber() != resnumber:
					if resnumber in structure.residue_dict:
						curr_residue = Residue(resname, resnumber+10000)
						structure.add_residue(curr_residue)	
					else:
						curr_residue = Residue(resname, resnumber)
						structure.add_residue(curr_residue)	
			except NameError:
				if resnumber in structure.residue_dict:
					curr_residue = Residue(resname, resnumber+10000)
					structure.add_residue(curr_residue)	
				else:
					curr_residue = Residue(resname, resnumber)
					structure.add_residue(curr_residue)	
			structure.residue_dict[resnumber].add_atom(curr_atom)
		elif record_type == "TER":
			curr_ter = Ter_record(resnumber)
			structure.add_ter(curr_ter)
		elif record_type == "CONECT":
			atom1, atom2 = line[7:20].split()
			curr_conect = Conect_record(atom1, atom2)
			structure.add_conect(curr_conect)
		elif record_type == "CRYST1":
			xbox = float(line[8:15])
			ybox = float(line[17:24])
			zbox = float(line[26:33])
			xangle = float(line[34:40])
			yangle = float(line[41:47])
			zangle = float(line[48:54])
			cryst1 = Cryst1_record(xbox, ybox, zbox, xangle, yangle, zangle)
			structure.add_cryst1(cryst1)
	return structure

def savepdb(struct, out):
	pdb = open(out, 'w')
	try:
		pdb.write(struct.other_dict['Cryst1'].formatted())
	except KeyError:
		pass
	for res in struct.residue_list:
		for atom in res.atom_list:
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
	pdb.close()

def sulpher(struct, pdb):
	sulfurs = []
	for res in struct.residue_list:
		if res.get_resname() == 'CYX':
			sulfurs.append(res.atom_dict['SG'])
		if res.get_resname() in [ 'OEH', 'OEC', 'OET', 'I7X' ] :
			for res2 in struct.residue_list:
				if res2.get_resname() in ['MKC','MKT','NL7', 'I4X']:
					curr_conect = Conect_record(res.atom_dict['CT'].get_number(), res2.atom_dict['CE'].get_number())
					struct.add_conect(curr_conect)
		if res.get_resname() == 'NLS':
			res2 = res.get_resnumber() + 4
			if struct.residue_dict[res2].get_resname() in ['NLS']:
				curr_conect = Conect_record(res.atom_dict['CE'].get_number(), struct.residue_dict[res2].atom_dict['CE'].get_number())
	                        struct.add_conect(curr_conect)
	for sul1 in sulfurs[:]:
		s1 = {}
		try:
			sulfurs.remove(sul1)
			sulph = sulfurs
			for sul2 in sulph[:]:
				key, value = sul2, sul1.calc_distance(sul2)
				s1[key] = float(value)
			prox = min(s1, key=s1.get)
			curr_conect = Conect_record(sul1.get_number(), prox.get_number())
			struct.add_conect(curr_conect)
			sulfurs.remove(prox)
		except:
			pass
	pdb = open(pdb, 'w')
	try:
		pdb.write(struct.other_dict['Cryst1'].formatted())
	except KeyError:
		pass
	for res in struct.residue_list:
		for atom in res.atom_list:
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
	pdb.close()

def stapler(struct, kind):
	if kind == 'OEH':
		for res in struct.residue_list:
			if res.get_resname() == 'NX1':
				nx1 = res
		nx1.set_resname('OEH')
		for atom in nx1.atom_list[:]:
		        if atom.get_name() not in ['N', 'H', 'CA', 'CB2', 'HB21', 'HB22', 'HB23', 'CB1', 'HB12', 'HB13', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'CH', 'HH2', 'HH3', 'CT', 'HT', 'C', 'O']:
                                nx1.atom_list.remove(atom)
                for res in struct.residue_list:
                        if res.get_resname() == 'NX2':
                                nx2 = res
		nx2.set_resname('NL7')
		for atom in nx2.atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'CB2', 'HB21', 'HB22', 'HB23', 'CB1', 'HB12', 'HB13', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE', 'C', 'O']:
                                nx2.atom_list.remove(atom)
		        pdb = open('Mutated.pdb', 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
                        if atom.get_name() == 'PC':
                                atom.set_name('Cl-')
                                atom.set_resname('Cl-')
                        if atom.get_name() == 'PN':
                                atom.set_name('Na+')
                                atom.set_resname('Na+')
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
        pdb.close()
        ctrl = open('chop.in', 'w')
        ctrl.write("source leaprc.ff99SB+\n")
        ctrl.write("Mut = loadpdb Mutated.pdb\n")
        ctrl.write("savepdb Mut Mut_leap.pdb\n")
        ctrl.write("quit\n")
        ctrl.close()
        Leapy.run('chop.in')
        if struct.other_dict['Cryst1']:
                f = open('Mut_leap.pdb', 'r+')
                temp = f.read()
                f.close()
                f = open('Mut_leap.pdb', 'w+')
                f.write(struct.other_dict['Cryst1'].formatted())
                f.write(temp)
                f.close()
        s = readpdb('Mut_leap.pdb')
        sulpher(s, 'Mut_leap.pdb')
        os.remove('chop.in')


def chop(struct, resid, resnum):
	struct.residue_dict[resnum].set_resname(resid)
	if resid == 'ABU':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG1', 'HG2', 'HG3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid in ('ALA', 'DAA'):
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid in ('ALX'):
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'ANN':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE', 'HE',  'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'ARG':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'AMX':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'ND2', 'HD2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'ASN':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'ASP':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'CBA':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'CD2', 'HD21', 'HD21', 'CE', 'HE2', 'HE3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'CYS':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'HG', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'GLN':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'GLU':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'OE2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'GLY':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA2', 'HA3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HCY':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'HD', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HCX':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HIS':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'HE2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HR1':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ', 'CH', 'NT1', 'HT11', 'HT12', 'NT2', 'HT21', 'HT22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HR2':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'NH', 'HH', 'CT', 'NI1', 'HI11', 'HI12', 'NI2', 'HI21', 'HI22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HR3':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'CH', 'HH2', 'HH3', 'NT', 'HT', 'CI', 'NK1', 'HK11', 'HK12', 'NK2', 'HK21', 'HK22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HR4':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'CH', 'HH2', 'HH3', 'CT', 'HT2', 'HT3', 'NI', 'HI', 'CK', 'NL1', 'HL11', 'HL12', 'NL2', 'HL21', 'HL22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'HR5':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'CH', 'HH2', 'HH3', 'CT', 'HT2', 'HT3', 'CI', 'HI2', 'HI3', 'NK', 'HK', 'CL', 'NM1', 'HM11', 'HM12', 'NM2', 'HM21', 'HM22', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid in ['I4X']:
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'CB2', 'HB21', 'HB22', 'HB23', 'CB1', 'HB12', 'HB13', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
		struct.residue_dict[resnum + 4].set_resname(resid)
		for atom in struct.residue_dict[resnum + 4].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'CB2', 'HB21', 'HB22', 'HB23', 'CB1', 'HB12', 'HB13', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE', 'C', 'O']:
				struct.residue_dict[resnum + 4].atom_list.remove(atom)	
	if resid == 'I7X':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'CB2', 'HB21', 'HB22', 'HB23', 'CB1', 'HB12', 'HB13', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'CH', 'HH2', 'HH3', 'CT', 'HT', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'ILE':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'HG13', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'LEU':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'LYS':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'MCX':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'CD', 'HD2', 'HD3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'MCY':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'CD', 'HD1', 'HD2', 'HD3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'MET':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'NER':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'CH', 'HH1', 'HH2', 'HH3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'NKI':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ2', 'HZ3', 'CH', 'HH2', 'HH3', 'CT', 'HT1', 'HT2', 'HT3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'NLE':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'NMA':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'CZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'NVA':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD1', 'HD2', 'HD3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'PHE':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'HZ', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'PRO':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'QUA':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'SER':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'OG', 'HG', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'THR':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'OG1', 'HG1', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'TRP':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'TYR':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'OH', 'HH', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'TZ1':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'NG', 'ND1', 'NE1', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'TZ2':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'ND', 'NE1', 'NZ1', 'CZ2', 'HZ2', 'CE2', 'HE2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'TZ3':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'NZ1', 'NH1', 'CH2', 'HH2', 'CZ2', 'HZ2', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	if resid == 'VAL':
		for atom in struct.residue_dict[resnum].atom_list[:]:
			if atom.get_name() not in ['N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O']:
				struct.residue_dict[resnum].atom_list.remove(atom)	
	pdb = open('Mutated.pdb', 'w')
	try:
		pdb.write(struct.other_dict['Cryst1'].formatted())
	except KeyError:
		pass
	for res in struct.residue_list:
		for atom in res.atom_list:
			if atom.get_name() == 'PC':
				atom.set_name('Cl-')
				atom.set_resname('Cl-')
			if atom.get_name() == 'PN':
				atom.set_name('Na+')
				atom.set_resname('Na+')
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
	pdb.close()
	ctrl = open('chop.in', 'w')
        ctrl.write("source leaprc.ff99SB+\n")
        ctrl.write("source leaprc.lipid14\n")
	ctrl.write("Mut = loadpdb Mutated.pdb\n")
	ctrl.write("savepdb Mut Mut_leap.pdb\n")
	ctrl.write("quit\n")
	ctrl.close()
	Leapy.run('chop.in')
	if struct.other_dict['Cryst1']:
		f = open('Mut_leap.pdb', 'r+')
		temp = f.read()
		f.close()
		f = open('Mut_leap.pdb', 'w+')
		f.write(struct.other_dict['Cryst1'].formatted())
		f.write(temp)
		f.close()
	s = readpdb('Mut_leap.pdb')
	sulpher(s, 'Mut_leap.pdb')
	os.remove('chop.in')
