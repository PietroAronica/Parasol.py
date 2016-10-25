import numpy

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

	def halfway_between(self, name, atom1, atom2):
		x = (atom2.get_coord()[0] + atom1.get_coord()[0])/2
		y = (atom2.get_coord()[1] + atom1.get_coord()[1])/2
		z = (atom2.get_coord()[2] + atom1.get_coord()[2])/2
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
					curr.residue = Residue(resname, resnumber)
					structure.add_residue(curr_residue)	
			except NameError:
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

