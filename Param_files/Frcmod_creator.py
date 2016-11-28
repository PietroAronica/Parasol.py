import fileinput

def make(file):
	fmod=open(file, 'w+')
	fmod.write("#Forcemod File" + '\n')
	fmod.write("MASS" + '\n')
	fmod.write('\n')
	fmod.write("BOND" + '\n')
	fmod.write('\n')
	fmod.write("ANGLE" + '\n')
	fmod.write('\n')
	fmod.write("DIHEDRAL" + '\n')
	fmod.write('\n')
	fmod.write("IMPROPER" + '\n')
	fmod.write('\n')
	fmod.write("NONBON" + '\n')
	fmod.write('\n')
	fmod.close()

def make_hyb():
	hyb=open('Hyb.dat', 'w+')
	hyb.write("addAtomTypes {" + '\n')
	hyb.write("}" + '\n')
	hyb.close()

def TYPE_insert(typ, ele, hyb):
	proc= False
	for line in fileinput.input('Hyb.dat', inplace=1):
		if line.startswith('addAtomTypes'):
			proc = True
		else:
			if proc:
				print '{{ \"{}\"  \"{}\" \"{}\" }}'.format(typ, ele, hyb)
			proc = False
		print line,

def MASS_insert(file, atom, mass, pola):
	proc= False
	for line in fileinput.input(file, inplace=1):
		if line.startswith('MASS'):
			proc = True
		else:
			if proc:
				print atom, mass, pola
			proc = False
		print line,

def BOND_insert(file, bond, fcst, leng):
	proc= False
	for line in fileinput.input(file, inplace=1):
		if line.startswith('BOND'):
			proc = True
		else:
			if proc:
				print bond, fcst, leng
			proc = False
		print line,

def ANGLE_insert(file, angl, fcst, widt):
	proc= False
	for line in fileinput.input(file, inplace=1):
		if line.startswith('ANGLE'):
			proc = True
		else:
			if proc:
				print angl, fcst, widt
			proc = False
		print line,

def DIHEDRAL_insert(file, dihe, divi, fcst, angl, peri):
	proc= False
	for line in fileinput.input(file, inplace=1):
		if line.startswith('DIHEDRAL'):
			proc = True
		else:
			if proc:
				print '{} {:.0f} {} {} {:.0f}.'.format(dihe, divi, fcst, angl, peri)
			proc = False
		print line,

def IMPROPER_insert(file, impr, fcst, angl, peri):
	proc= False
	for line in fileinput.input(file, inplace=1):
		if line.startswith('IMPROPER'):
			proc = True
		else:
			if proc:
				print '{}    {} {} {:.0f}.'.format(impr, fcst, angl, peri)
			proc = False
		print line,

def NONBON_insert(file, atom, vdwr, potw):
	proc= False
	for line in fileinput.input(file, inplace=1):
		if line.startswith('NONBON'):
			proc = True
		else:
			if proc:
				print atom, vdwr, potw
			proc = False
		print line,

